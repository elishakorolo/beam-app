import React, { useState, useRef, useEffect } from 'react';
import { AreaChart, Area, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer } from 'recharts';

// ============================================================================
// FRAME ANALYSIS ENGINE
// ============================================================================

class FrameAnalyzer {
  constructor(nodes, members, canSway) {
    this.nodes = JSON.parse(JSON.stringify(nodes));
    this.members = JSON.parse(JSON.stringify(members));
    this.memberResults = new Map();
    this.canSway = canSway;
  }

  getMemberGeometry(m) {
    const nodeA = this.nodes.find(n => n.id === m.start);
    const nodeB = this.nodes.find(n => n.id === m.end);
    const dx = nodeB.x - nodeA.x;
    const dy = nodeB.y - nodeA.y;
    const L = Math.sqrt(dx * dx + dy * dy);
    const angle = Math.atan2(dy, dx);
    const cos = Math.cos(angle);
    const sin = Math.sin(angle);
    return { nodeA, nodeB, L, angle, cos, sin, isVert: Math.abs(cos) < 0.01, isHoriz: Math.abs(sin) < 0.01 };
  }

  calcFEM(m, geom) {
    let femAB = 0, femBA = 0;
    const L = geom.L;

    for (const load of m.loads) {
      if (load.type === 'point' && load.pos !== undefined) {
        const a = load.pos;
        const b = L - a;
        femAB -= (load.mag * a * b * b) / (L * L);
        femBA += (load.mag * a * a * b) / (L * L);
      } else if (load.type === 'udl') {
        femAB -= (load.mag * L * L) / 12;
        femBA += (load.mag * L * L) / 12;
      } else if (load.type === 'vdl' && load.magEnd !== undefined) {
        const w1 = load.mag;
        const w2 = load.magEnd;
        const wRect = Math.min(w1, w2);
        const wTri = Math.abs(w2 - w1);
        femAB -= (wRect * L * L) / 12;
        femBA += (wRect * L * L) / 12;
        if (w2 > w1) {
          femAB -= (wTri * L * L) / 20;
          femBA += (wTri * L * L) / 30;
        } else {
          femAB -= (wTri * L * L) / 30;
          femBA += (wTri * L * L) / 20;
        }
      }
    }
    return [femAB, femBA];
  }

  solve() {
    const geoms = new Map();
    for (const m of this.members) {
      const geom = this.getMemberGeometry(m);
      const [femAB, femBA] = this.calcFEM(m, geom);
      geoms.set(m.id, { ...geom, femAB, femBA });
    }

    const rotDofs = this.nodes.filter(n => n.support !== 'fixed');
    const dofMap = new Map(rotDofs.map((n, i) => [n.id, i]));
    const nDof = rotDofs.length + (this.canSway ? 1 : 0);

    if (nDof === 0) throw new Error('Fully constrained');

    const K = Array(nDof).fill(0).map(() => Array(nDof).fill(0));
    const F = Array(nDof).fill(0);

    for (const m of this.members) {
      const g = geoms.get(m.id);
      const I_eff = m.I * (m.I_mult || 1);
      const k = 2 * m.E * I_eff / g.L;
      const iA = dofMap.get(m.start);
      const iB = dofMap.get(m.end);

      if (iA !== undefined) {
        K[iA][iA] += 2 * k;
        F[iA] -= g.femAB;
        if (iB !== undefined) K[iA][iB] += k;
        if (this.canSway && g.isVert) K[iA][nDof - 1] -= 3 * k / g.L;
      }

      if (iB !== undefined) {
        K[iB][iB] += 2 * k;
        F[iB] -= g.femBA;
        if (iA !== undefined) K[iB][iA] += k;
        if (this.canSway && g.isVert) K[iB][nDof - 1] -= 3 * k / g.L;
      }
    }

    if (this.canSway) {
      for (const m of this.members) {
        const g = geoms.get(m.id);
        if (!g.isVert) continue;
        const I_eff = m.I * (m.I_mult || 1);
        const k = 2 * m.E * I_eff / g.L;
        F[nDof - 1] -= (g.femAB + g.femBA) / g.L;
        const iA = dofMap.get(m.start);
        const iB = dofMap.get(m.end);
        if (iA !== undefined) K[nDof - 1][iA] += 3 * k / g.L;
        if (iB !== undefined) K[nDof - 1][iB] += 3 * k / g.L;
        K[nDof - 1][nDof - 1] -= 6 * k / (g.L * g.L);
      }
    }

    const x = this.gaussElim(K, F);
    rotDofs.forEach((n, i) => n.theta = x[i]);
    const sway = this.canSway ? x[nDof - 1] : 0;

    for (const m of this.members) {
      const g = geoms.get(m.id);
      const I_eff = m.I * (m.I_mult || 1);
      const k = 2 * m.E * I_eff / g.L;
      const tA = g.nodeA.theta;
      const tB = g.nodeB.theta;
      const psi = (this.canSway && g.isVert) ? sway / g.L : 0;
      const mAB = g.femAB + k * (2 * tA + tB - 3 * psi);
      const mBA = g.femBA + k * (tA + 2 * tB - 3 * psi);

      let totTrans = 0, momB = 0;
      for (const load of m.loads) {
        if (load.type === 'point' && load.pos !== undefined) {
          totTrans += load.mag;
          momB += load.mag * (g.L - load.pos);
        } else if (load.type === 'udl') {
          const tot = load.mag * g.L;
          totTrans += tot;
          momB += tot * g.L / 2;
        } else if (load.type === 'vdl' && load.magEnd !== undefined) {
          const avg = (load.mag + load.magEnd) / 2;
          const tot = avg * g.L;
          totTrans += tot;
          let centroid = g.L / 2;
          if (Math.abs(load.magEnd - load.mag) > 1e-6) {
            if (load.magEnd > load.mag) {
              centroid = g.L * (load.mag + 2 * load.magEnd) / (3 * (load.mag + load.magEnd));
            } else {
              centroid = g.L * (2 * load.mag + load.magEnd) / (3 * (load.mag + load.magEnd));
            }
          }
          momB += tot * (g.L - centroid);
        }
      }
      const vA = (momB - mAB - mBA) / g.L;
      const vB = -(totTrans - vA);

      this.memberResults.set(m.id, {
        L: g.L, angle: g.angle, cos: g.cos, sin: g.sin,
        isVert: g.isVert, isHoriz: g.isHoriz,
        femAB: g.femAB, femBA: g.femBA,
        mAB, mBA, vA, vB
      });
    }

    for (const node of this.nodes) {
      if (node.support === 'free') continue;
      let fx = 0, fy = 0, m = 0;
      for (const mem of this.members) {
        const res = this.memberResults.get(mem.id);
        if (mem.start === node.id) {
          fx += res.vA * res.sin;
          fy += -res.vA * res.cos;
          m += res.mAB;
        } else if (mem.end === node.id) {
          fx += -res.vB * res.sin;
          fy += res.vB * res.cos;
          m += res.mBA;
        }
      }
      node.rx = -fx;
      node.ry = -fy;
      if (node.support === 'fixed') node.rm = -m;
    }

    return sway;
  }

  gaussElim(A, b) {
    const n = A.length;
    const aug = A.map((row, i) => [...row, b[i]]);

    for (let i = 0; i < n; i++) {
      let max = i;
      for (let k = i + 1; k < n; k++) {
        if (Math.abs(aug[k][i]) > Math.abs(aug[max][i])) max = k;
      }
      [aug[i], aug[max]] = [aug[max], aug[i]];
      if (Math.abs(aug[i][i]) < 1e-10) throw new Error('Singular matrix');
      for (let k = i + 1; k < n; k++) {
        const f = aug[k][i] / aug[i][i];
        for (let j = i; j <= n; j++) aug[k][j] -= f * aug[i][j];
      }
    }

    const x = Array(n).fill(0);
    for (let i = n - 1; i >= 0; i--) {
      x[i] = aug[i][n];
      for (let j = i + 1; j < n; j++) x[i] -= aug[i][j] * x[j];
      x[i] /= aug[i][i];
    }
    return x;
  }

  getDiagram(memId, type, pts = 100) {
    const m = this.members.find(mm => mm.id === memId);
    const res = this.memberResults.get(memId);
    const data = [];

    for (let i = 0; i <= pts; i++) {
      const s = (res.L * i) / pts;
      let val = 0;

      if (type === 'bmd') {
        val = res.mAB + res.vA * s;
        
        for (const load of m.loads) {
          if (load.type === 'point' && load.pos !== undefined && s >= load.pos) {
            val -= load.mag * (s - load.pos);
          } else if (load.type === 'udl') {
            val -= 0.5 * load.mag * s * s;
          } else if (load.type === 'vdl' && load.magEnd !== undefined) {
            const w1 = load.mag;
            const w2 = load.magEnd;
            const w_s = w1 + (w2 - w1) * s / res.L;
            const avg_load = (w1 + w_s) / 2;
            val -= 0.5 * avg_load * s * s;
          }
        }
      } else {
        val = res.vA;
        
        for (const load of m.loads) {
          if (load.type === 'point' && load.pos !== undefined && s >= load.pos) {
            val -= load.mag;
          } else if (load.type === 'udl') {
            val -= load.mag * s;
          } else if (load.type === 'vdl' && load.magEnd !== undefined) {
            const w1 = load.mag;
            const w2 = load.magEnd;
            const w_s = w1 + (w2 - w1) * s / res.L;
            const avg_load = (w1 + w_s) / 2;
            val -= avg_load * s;
          }
        }
      }
      
      data.push({ s, val });
    }
    return data;
  }
}

// ============================================================================
// DESIGN FUNCTIONS
// ============================================================================

function designBeam(M, fck, fy, b, cover) {
  const bMm = b * 1000;
  const fcd = 0.67 * fck / 1.5;
  const fyd = fy / 1.15;
  const xuLim = 0.48;
  const MuLim = 0.36 * xuLim * (1 - 0.42 * xuLim) * fcd;
  const Md = Math.abs(M) * 1e6;
  const dReq = Math.sqrt(Md / (MuLim * bMm));
  let D = Math.ceil((dReq + cover + 16) / 25) * 25;
  const d = D - cover - 16;
  const MuLimAct = MuLim * bMm * d * d;
  let Ast = 0, type = '';
  if (Md <= MuLimAct) {
    Ast = 0.5 * fcd * bMm * d * (1 - Math.sqrt(1 - 4.6 * Md / (fcd * bMm * d * d))) / fyd;
    type = 'Singly';
  } else {
    Ast = 0.5 * fcd * bMm * d / fyd + (Md - MuLimAct) / (fyd * 0.9 * d);
    type = 'Doubly';
  }
  Ast = Math.max(Ast, 0.85 * bMm * d / fy);
  return { D, d, Ast, type };
}

function designColumn(P, M, fck, fy, size, cover) {
  const fcd = 0.67 * fck / 1.5;
  const fyd = fy / 1.15;
  const PuMax = 0.4 * fck * size * size / 1000;
  let Asc = P > 0.1 * PuMax 
    ? Math.max(0.008 * size * size, (P * 1000 - 0.4 * fcd * size * size) / fyd)
    : 0.008 * size * size;
  Asc = Math.max(Asc, 0.008 * size * size);
  Asc = Math.min(Asc, 0.04 * size * size);
  return { size, Asc, pct: (Asc / (size * size)) * 100, PuMax, util: PuMax > 0 ? (P / PuMax) * 100 : 0 };
}

// ============================================================================
// SVG FRAME RENDERER - STRUCTURE ONLY
// ============================================================================

function FrameView({ nodes, members, results, mode }) {
  const svgRef = useRef(null);
  const [viewBox, setViewBox] = useState('0 0 800 600');

  useEffect(() => {
    if (nodes.length === 0) return;

    const xs = nodes.map(n => n.x);
    const ys = nodes.map(n => n.y);
    const minX = Math.min(...xs);
    const maxX = Math.max(...xs);
    const minY = Math.min(...ys);
    const maxY = Math.max(...ys);

    const w = maxX - minX || 10;
    const h = maxY - minY || 10;
    const pad = Math.max(w, h) * 0.2;

    setViewBox(`${minX - pad} ${minY - pad} ${w + 2 * pad} ${h + 2 * pad}`);
  }, [nodes]);

  const renderSupport = (n) => {
    const size = 0.3;
    if (n.support === 'pinned') {
      return (
        <g key={`sup-${n.id}`}>
          <polygon
            points={`${n.x},${n.y} ${n.x - size},${n.y + size * 1.5} ${n.x + size},${n.y + size * 1.5}`}
            fill="#60a5fa"
            stroke="#1e40af"
            strokeWidth="0.05"
          />
          <line x1={n.x - size * 1.2} y1={n.y + size * 1.5} x2={n.x + size * 1.2} y2={n.y + size * 1.5} stroke="#1e40af" strokeWidth="0.1" />
        </g>
      );
    } else if (n.support === 'fixed') {
      return (
        <g key={`sup-${n.id}`}>
          <rect x={n.x - size} y={n.y - size} width={size * 2} height={size * 2} fill="#f59e0b" stroke="#b45309" strokeWidth="0.05" />
          {[...Array(4)].map((_, i) => (
            <line
              key={i}
              x1={n.x - size}
              y1={n.y - size + i * size / 2}
              x2={n.x + size}
              y2={n.y - size + i * size / 2}
              stroke="#78350f"
              strokeWidth="0.02"
            />
          ))}
        </g>
      );
    } else if (n.support === 'roller') {
      return (
        <g key={`sup-${n.id}`}>
          <polygon
            points={`${n.x},${n.y} ${n.x - size},${n.y + size * 1.2} ${n.x + size},${n.y + size * 1.2}`}
            fill="#10b981"
            stroke="#065f46"
            strokeWidth="0.05"
          />
          <circle cx={n.x - size * 0.5} cy={n.y + size * 1.5} r={size * 0.2} fill="white" stroke="#065f46" strokeWidth="0.03" />
          <circle cx={n.x + size * 0.5} cy={n.y + size * 1.5} r={size * 0.2} fill="white" stroke="#065f46" strokeWidth="0.03" />
        </g>
      );
    }
    return null;
  };

  return (
    <svg ref={svgRef} viewBox={viewBox} style={{ width: '100%', height: '100%', transform: 'scaleY(-1)' }}>
      {/* Title */}
      <g transform={`scale(1, -1)`}>
        <rect
          x={nodes.length > 0 ? Math.min(...nodes.map(n => n.x)) - 0.5 : 0}
          y={-(nodes.length > 0 ? Math.max(...nodes.map(n => n.y)) + 1.2 : 5)}
          width="3"
          height="0.4"
          fill="rgba(255, 255, 255, 0.95)"
          stroke="#e5e7eb"
          strokeWidth="0.03"
          rx="0.1"
        />
        <text
          x={nodes.length > 0 ? Math.min(...nodes.map(n => n.x)) - 0.3 : 0.2}
          y={-(nodes.length > 0 ? Math.max(...nodes.map(n => n.y)) + 0.9 : 4.8)}
          fontSize="0.2"
          fill="#1f2937"
          fontWeight="bold"
          fontFamily="sans-serif"
        >
          Frame Structure
        </text>
      </g>
      
      {/* Members */}
      {members.map(m => {
        const nA = nodes.find(n => n.id === m.start);
        const nB = nodes.find(n => n.id === m.end);
        return (
          <line
            key={`mem-${m.id}`}
            x1={nA.x}
            y1={nA.y}
            x2={nB.x}
            y2={nB.y}
            stroke="#94a3b8"
            strokeWidth="0.15"
            strokeLinecap="round"
          />
        );
      })}

      {/* Member Labels */}
      {members.map(m => {
        const nA = nodes.find(n => n.id === m.start);
        const nB = nodes.find(n => n.id === m.end);
        const midX = (nA.x + nB.x) / 2;
        const midY = (nA.y + nB.y) / 2;
        return (
          <g key={`mem-label-${m.id}`} transform={`scale(1, -1)`}>
            <circle cx={midX} cy={-midY} r="0.2" fill="white" stroke="#3b82f6" strokeWidth="0.04" />
            <text
              x={midX}
              y={-midY + 0.07}
              fontSize="0.15"
              fill="#3b82f6"
              fontWeight="bold"
              textAnchor="middle"
              fontFamily="monospace"
            >
              {m.id + 1}
            </text>
          </g>
        );
      })}

      {/* Supports */}
      {nodes.map(renderSupport)}

      {/* Nodes */}
      {nodes.map(n => (
        <g key={`node-${n.id}`}>
          <circle cx={n.x} cy={n.y} r="0.15" fill="#f59e0b" stroke="#b45309" strokeWidth="0.04" />
          <g transform={`scale(1, -1)`}>
            <text
              x={n.x + 0.25}
              y={-n.y + 0.05}
              fontSize="0.18"
              fill="#1f2937"
              fontWeight="bold"
              fontFamily="sans-serif"
            >
              {n.id + 1}
            </text>
          </g>
        </g>
      ))}
    </svg>
  );
}

// ============================================================================
// STYLES
// ============================================================================

const styles = {
  container: {
    minHeight: '100vh',
    background: 'linear-gradient(to bottom right, #f9fafb, #e5e7eb)',
    color: '#111827',
    fontFamily: 'system-ui, -apple-system, sans-serif'
  },
  header: {
    background: 'rgba(255, 255, 255, 0.8)',
    backdropFilter: 'blur(20px)',
    borderBottom: '1px solid #e5e7eb',
    position: 'sticky',
    top: 0,
    zIndex: 1000,
    boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)'
  },
  title: {
    fontSize: '1.5rem',
    fontWeight: '300',
    color: '#111827',
    letterSpacing: '-0.025em'
  },
  subtitle: {
    fontSize: '0.75rem',
    color: '#6b7280',
    marginTop: '0.25rem'
  },
  analyzeBtn: {
    padding: '0.625rem 1.5rem',
    background: '#3b82f6',
    color: 'white',
    fontSize: '0.875rem',
    fontWeight: '500',
    borderRadius: '0.5rem',
    border: 'none',
    cursor: 'pointer',
    boxShadow: '0 4px 6px -1px rgba(59, 130, 246, 0.3)',
    transition: 'all 0.3s'
  },
  sidebar: {
    width: '20rem',
    background: 'rgba(255, 255, 255, 0.8)',
    backdropFilter: 'blur(20px)',
    borderRight: '1px solid #e5e7eb',
    padding: '1.5rem',
    overflowY: 'auto',
    height: 'calc(100vh - 80px)',
    position: 'sticky',
    top: '80px'
  },
  sectionTitle: {
    fontSize: '0.75rem',
    fontWeight: '600',
    color: '#374151',
    marginBottom: '0.75rem',
    textTransform: 'uppercase',
    letterSpacing: '0.05em'
  },
  input: {
    width: '100%',
    padding: '0.5rem 0.75rem',
    background: '#f9fafb',
    border: '1px solid #d1d5db',
    borderRadius: '0.5rem',
    color: '#111827',
    fontSize: '0.875rem',
    outline: 'none',
    transition: 'border-color 0.2s'
  },
  label: {
    fontSize: '0.75rem',
    color: '#6b7280',
    marginBottom: '0.25rem',
    display: 'block'
  },
  select: {
    width: '100%',
    padding: '0.5rem',
    background: '#f9fafb',
    border: '1px solid #d1d5db',
    borderRadius: '0.5rem',
    color: '#111827',
    fontSize: '0.875rem'
  },
  card: {
    padding: '0.75rem',
    background: 'white',
    border: '1px solid #e5e7eb',
    borderRadius: '0.75rem',
    marginBottom: '0.75rem',
    boxShadow: '0 1px 2px 0 rgba(0, 0, 0, 0.05)'
  },
  button: {
    padding: '0.375rem 0.75rem',
    background: '#dbeafe',
    color: '#1e40af',
    border: '1px solid #bfdbfe',
    borderRadius: '0.5rem',
    cursor: 'pointer',
    fontSize: '0.75rem',
    fontWeight: '500'
  },
  tab: {
    flex: 1,
    padding: '0.625rem',
    fontSize: '0.75rem',
    fontWeight: '500',
    border: 'none',
    background: 'transparent',
    color: '#6b7280',
    cursor: 'pointer',
    transition: 'all 0.2s',
    borderRadius: '0.5rem'
  },
  tabActive: {
    background: '#3b82f6',
    color: 'white',
    boxShadow: '0 1px 3px 0 rgba(59, 130, 246, 0.3)'
  },
  viewBtn: {
    padding: '0.5rem 1.5rem',
    fontSize: '0.875rem',
    fontWeight: '500',
    borderRadius: '0.5rem',
    border: 'none',
    cursor: 'pointer',
    transition: 'all 0.3s'
  },
  viewBtnActive: {
    background: '#3b82f6',
    color: 'white',
    boxShadow: '0 4px 6px -1px rgba(59, 130, 246, 0.3)'
  },
  viewBtnInactive: {
    background: 'white',
    color: '#6b7280',
    border: '1px solid #e5e7eb'
  },
  resultCard: {
    padding: '1.5rem',
    background: 'white',
    border: '1px solid #e5e7eb',
    borderRadius: '1rem',
    boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)'
  },
  gradientCard: {
    padding: '1rem',
    borderRadius: '1rem',
    boxShadow: '0 1px 3px 0 rgba(0, 0, 0, 0.1)'
  }
};

// ============================================================================
// MAIN COMPONENT
// ============================================================================

export default function App() {
  const [E, setE] = useState(200);
  const [nodes, setNodes] = useState([
    { id: 0, x: 0, y: 0, support: 'fixed', theta: 0, rx: 0, ry: 0, rm: 0 },
    { id: 1, x: 0, y: 4, support: 'free', theta: 0, rx: 0, ry: 0, rm: 0 },
    { id: 2, x: 5, y: 4, support: 'free', theta: 0, rx: 0, ry: 0, rm: 0 },
    { id: 3, x: 5, y: 0, support: 'pinned', theta: 0, rx: 0, ry: 0, rm: 0 }
  ]);
  const [members, setMembers] = useState([
    { id: 0, start: 0, end: 1, E: 200e9, I: 1e-6, I_mult: 1, loads: [] },
    { id: 1, start: 1, end: 2, E: 200e9, I: 1e-6, I_mult: 1, loads: [{ type: 'udl', mag: 10 }] },
    { id: 2, start: 2, end: 3, E: 200e9, I: 1e-6, I_mult: 1, loads: [] }
  ]);
  const [sway, setSway] = useState(false);
  const [fck, setFck] = useState(25);
  const [fy, setFy] = useState(415);
  const [beamW, setBeamW] = useState(0.3);
  const [colSize, setColSize] = useState(300);
  const [cover, setCover] = useState(40);

  const [results, setResults] = useState(null);
  const [view, setView] = useState('structure');
  const [activePanel, setActivePanel] = useState('nodes');

  const addNode = () => {
    const newId = nodes.length;
    setNodes([...nodes, { id: newId, x: 0, y: 0, support: 'free', theta: 0, rx: 0, ry: 0, rm: 0 }]);
  };

  const updateNode = (index, field, value) => {
    const newNodes = [...nodes];
    newNodes[index][field] = value;
    setNodes(newNodes);
  };

  const deleteNode = (index) => {
    if (nodes.length <= 2) {
      alert('Need at least 2 nodes!');
      return;
    }
    const newNodes = nodes.filter((_, i) => i !== index);
    newNodes.forEach((n, i) => n.id = i);
    setNodes(newNodes);
  };

  const addMember = () => {
    if (nodes.length < 2) return;
    const newId = members.length;
    setMembers([...members, { id: newId, start: 0, end: 1, E: E * 1e9, I: 1e-6, I_mult: 1, loads: [] }]);
  };

  const updateMember = (index, field, value) => {
    const newMembers = [...members];
    newMembers[index][field] = value;
    setMembers(newMembers);
  };

  const deleteMember = (index) => {
    const newMembers = members.filter((_, i) => i !== index);
    newMembers.forEach((m, i) => m.id = i);
    setMembers(newMembers);
  };

  const addLoad = (mIdx) => {
    const newMembers = [...members];
    newMembers[mIdx].loads.push({ type: 'udl', mag: 10 });
    setMembers(newMembers);
  };

  const updateLoad = (mIdx, lIdx, field, value) => {
    const newMembers = [...members];
    newMembers[mIdx].loads[lIdx][field] = value;
    setMembers(newMembers);
  };

  const deleteLoad = (mIdx, lIdx) => {
    const newMembers = [...members];
    newMembers[mIdx].loads = newMembers[mIdx].loads.filter((_, i) => i !== lIdx);
    setMembers(newMembers);
  };

  const analyze = () => {
    try {
      const solver = new FrameAnalyzer(nodes, members, sway);
      const swayDelta = solver.solve();

      const maxMB = Math.max(
        ...Array.from(solver.memberResults.values())
          .filter(r => r.isHoriz)
          .flatMap(r => [Math.abs(r.mAB), Math.abs(r.mBA)]),
        0
      );
      
      const maxMC = Math.max(
        ...Array.from(solver.memberResults.values())
          .filter(r => r.isVert)
          .flatMap(r => [Math.abs(r.mAB), Math.abs(r.mBA)]),
        0
      );

      const beamDes = maxMB > 0 ? designBeam(maxMB, fck, fy, beamW, cover) : null;
      const colDes = maxMC > 0 ? designColumn(0, maxMC, fck, fy, colSize, cover) : null;

      setResults({
        analyzer: solver,
        nodes: solver.nodes,
        memberResults: solver.memberResults,
        swayDelta,
        maxMB,
        maxMC,
        beamDes,
        colDes
      });
    } catch (e) {
      alert(`Error: ${e.message}`);
    }
  };

  return (
    <div style={styles.container}>
      {/* Header */}
      <header style={styles.header}>
        <div style={{ padding: '1rem 1.5rem', display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
          <div>
            <h1 style={styles.title}>Frame Analyzer</h1>
            <p style={styles.subtitle}>Slope-Deflection Analysis & Design - Clear Diagrams</p>
          </div>
          <button
            onClick={analyze}
            style={styles.analyzeBtn}
            onMouseEnter={(e) => e.target.style.transform = 'scale(1.05)'}
            onMouseLeave={(e) => e.target.style.transform = 'scale(1)'}
          >
            Analyze
          </button>
        </div>
      </header>

      <div style={{ display: 'flex' }}>
        {/* Left Sidebar */}
        <aside style={styles.sidebar}>
          {/* Material */}
          <section style={{ marginBottom: '1.5rem' }}>
            <h3 style={styles.sectionTitle}>Material</h3>
            <div>
              <label style={styles.label}>E (GPa)</label>
              <input
                type="number"
                value={E}
                onChange={e => setE(parseFloat(e.target.value) || 200)}
                style={styles.input}
              />
            </div>
          </section>

          {/* Tabs */}
          <div style={{ display: 'flex', gap: '0.25rem', background: '#f3f4f6', padding: '0.25rem', borderRadius: '0.5rem', marginBottom: '1.5rem' }}>
            {['nodes', 'members', 'loads'].map(panel => (
              <button
                key={panel}
                onClick={() => setActivePanel(panel)}
                style={{
                  ...styles.tab,
                  ...(activePanel === panel ? styles.tabActive : {})
                }}
              >
                {panel.charAt(0).toUpperCase() + panel.slice(1)}
              </button>
            ))}
          </div>

          {/* Nodes Panel */}
          {activePanel === 'nodes' && (
            <section>
              <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '1rem' }}>
                <h3 style={styles.sectionTitle}>Nodes</h3>
                <button onClick={addNode} style={styles.button}>+ Add</button>
              </div>
              <div style={{ display: 'flex', flexDirection: 'column', gap: '0.5rem', maxHeight: '400px', overflowY: 'auto' }}>
                {nodes.map((n, i) => (
                  <div key={n.id} style={styles.card}>
                    <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '0.5rem' }}>
                      <span style={{ fontSize: '0.75rem', color: '#3b82f6', fontWeight: '600' }}>Node {n.id + 1}</span>
                      {nodes.length > 2 && (
                        <button onClick={() => deleteNode(i)} style={{ ...styles.button, padding: '0.125rem 0.5rem', background: 'transparent', color: '#ef4444', border: 'none' }}>✕</button>
                      )}
                    </div>
                    <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '0.5rem', marginBottom: '0.5rem' }}>
                      <div>
                        <label style={styles.label}>X (m)</label>
                        <input
                          type="number"
                          value={n.x}
                          onChange={e => updateNode(i, 'x', parseFloat(e.target.value) || 0)}
                          style={styles.input}
                        />
                      </div>
                      <div>
                        <label style={styles.label}>Y (m)</label>
                        <input
                          type="number"
                          value={n.y}
                          onChange={e => updateNode(i, 'y', parseFloat(e.target.value) || 0)}
                          style={styles.input}
                        />
                      </div>
                    </div>
                    <div>
                      <label style={styles.label}>Support</label>
                      <select
                        value={n.support}
                        onChange={e => updateNode(i, 'support', e.target.value)}
                        style={styles.select}
                      >
                        <option value="free">Free</option>
                        <option value="pinned">Pinned</option>
                        <option value="fixed">Fixed</option>
                        <option value="roller">Roller</option>
                      </select>
                    </div>
                  </div>
                ))}
              </div>
            </section>
          )}

          {/* Members Panel */}
          {activePanel === 'members' && (
            <section>
              <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '1rem' }}>
                <h3 style={styles.sectionTitle}>Members</h3>
                <button onClick={addMember} style={styles.button}>+ Add</button>
              </div>
              <div style={{ display: 'flex', flexDirection: 'column', gap: '0.5rem', maxHeight: '400px', overflowY: 'auto' }}>
                {members.map((m, i) => (
                  <div key={m.id} style={styles.card}>
                    <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '0.5rem' }}>
                      <span style={{ fontSize: '0.75rem', color: '#3b82f6', fontWeight: '600' }}>Member {m.id + 1}</span>
                      <button onClick={() => deleteMember(i)} style={{ ...styles.button, padding: '0.125rem 0.5rem', background: 'transparent', color: '#ef4444', border: 'none' }}>✕</button>
                    </div>
                    <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '0.5rem' }}>
                      <div>
                        <label style={styles.label}>Start</label>
                        <select
                          value={m.start}
                          onChange={e => updateMember(i, 'start', parseInt(e.target.value))}
                          style={styles.select}
                        >
                          {nodes.map(n => <option key={n.id} value={n.id}>Node {n.id + 1}</option>)}
                        </select>
                      </div>
                      <div>
                        <label style={styles.label}>End</label>
                        <select
                          value={m.end}
                          onChange={e => updateMember(i, 'end', parseInt(e.target.value))}
                          style={styles.select}
                        >
                          {nodes.map(n => <option key={n.id} value={n.id}>Node {n.id + 1}</option>)}
                        </select>
                      </div>
                      <div style={{ gridColumn: 'span 2' }}>
                        <label style={styles.label}>I Multiplier</label>
                        <input
                          type="number"
                          value={m.I_mult || 1}
                          onChange={e => updateMember(i, 'I_mult', parseFloat(e.target.value) || 1)}
                          style={styles.input}
                          step="0.1"
                        />
                      </div>
                    </div>
                  </div>
                ))}
              </div>
            </section>
          )}

          {/* Loads Panel */}
          {activePanel === 'loads' && (
            <section>
              <h3 style={styles.sectionTitle}>Loads</h3>
              <div style={{ display: 'flex', flexDirection: 'column', gap: '0.5rem', maxHeight: '400px', overflowY: 'auto' }}>
                {members.map((m, mIdx) => (
                  <div key={m.id} style={styles.card}>
                    <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '0.5rem' }}>
                      <span style={{ fontSize: '0.75rem', color: '#3b82f6', fontWeight: '600' }}>Member {m.id + 1}</span>
                      <button onClick={() => addLoad(mIdx)} style={styles.button}>+</button>
                    </div>
                    {m.loads.map((l, lIdx) => (
                      <div key={lIdx} style={{ padding: '0.5rem', background: '#f0fdf4', borderRadius: '0.5rem', marginBottom: '0.5rem', border: '1px solid #bbf7d0' }}>
                        <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '0.5rem' }}>
                          <span style={{ fontSize: '0.75rem', color: '#16a34a' }}>Load {lIdx + 1}</span>
                          <button onClick={() => deleteLoad(mIdx, lIdx)} style={{ ...styles.button, padding: '0.125rem 0.5rem', background: 'transparent', color: '#ef4444', border: 'none' }}>✕</button>
                        </div>
                        <select
                          value={l.type}
                          onChange={e => updateLoad(mIdx, lIdx, 'type', e.target.value)}
                          style={{ ...styles.select, marginBottom: '0.5rem' }}
                        >
                          <option value="point">Point</option>
                          <option value="udl">UDL</option>
                          <option value="vdl">VDL</option>
                        </select>
                        <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '0.5rem' }}>
                          <div>
                            <label style={styles.label}>Mag {l.type !== 'point' ? '(kN/m)' : '(kN)'}</label>
                            <input
                              type="number"
                              value={l.mag}
                              onChange={e => updateLoad(mIdx, lIdx, 'mag', parseFloat(e.target.value) || 0)}
                              style={styles.input}
                            />
                          </div>
                          {l.type === 'point' && (
                            <div>
                              <label style={styles.label}>Pos (m)</label>
                              <input
                                type="number"
                                value={l.pos || 0}
                                onChange={e => updateLoad(mIdx, lIdx, 'pos', parseFloat(e.target.value) || 0)}
                                style={styles.input}
                              />
                            </div>
                          )}
                          {l.type === 'vdl' && (
                            <div>
                              <label style={styles.label}>End Mag</label>
                              <input
                                type="number"
                                value={l.magEnd || 0}
                                onChange={e => updateLoad(mIdx, lIdx, 'magEnd', parseFloat(e.target.value) || 0)}
                                style={styles.input}
                              />
                            </div>
                          )}
                        </div>
                      </div>
                    ))}
                  </div>
                ))}
              </div>
            </section>
          )}

          {/* Options */}
          <section style={{ marginTop: '1.5rem' }}>
            <h3 style={styles.sectionTitle}>Options</h3>
            <label style={{ display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
              <input type="checkbox" checked={sway} onChange={e => setSway(e.target.checked)} />
              <span style={{ fontSize: '0.875rem' }}>Allow Sway</span>
            </label>
          </section>

          {/* Design */}
          <section style={{ marginTop: '1.5rem' }}>
            <h3 style={styles.sectionTitle}>Design</h3>
            <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '0.5rem' }}>
              <div>
                <label style={styles.label}>fck (MPa)</label>
                <input type="number" value={fck} onChange={e => setFck(parseFloat(e.target.value) || 25)} style={styles.input} />
              </div>
              <div>
                <label style={styles.label}>fy (MPa)</label>
                <input type="number" value={fy} onChange={e => setFy(parseFloat(e.target.value) || 415)} style={styles.input} />
              </div>
            </div>
          </section>
        </aside>

        {/* Main Content */}
        <main style={{ flex: 1, padding: '2rem' }}>
          {/* View Selector */}
          <div style={{ display: 'flex', gap: '0.5rem', marginBottom: '1.5rem' }}>
            {['structure', 'bmd', 'sfd'].map(v => (
              <button
                key={v}
                onClick={() => setView(v)}
                style={{
                  ...styles.viewBtn,
                  ...(view === v ? styles.viewBtnActive : styles.viewBtnInactive)
                }}
              >
                {v.toUpperCase()}
              </button>
            ))}
          </div>

          {/* Frame Structure View - Always Show */}
          <div style={{ height: '400px', marginBottom: '1.5rem', background: 'white', borderRadius: '1rem', padding: '1rem', border: '1px solid #e5e7eb' }}>
            <FrameView nodes={nodes} members={members} results={results?.analyzer} mode="structure" />
          </div>

          {/* Member Diagrams - Show when BMD or SFD selected */}
          {results && view !== 'structure' && (
            <>
              {/* Diagram Summary */}
              <div style={{ ...styles.resultCard, marginBottom: '1.5rem', background: 'linear-gradient(135deg, #dbeafe, #f0f9ff)' }}>
                <h3 style={{ ...styles.sectionTitle, marginBottom: '1rem' }}>
                  {view === 'bmd' ? '📊 Bending Moment Diagrams' : '📊 Shear Force Diagrams'}
                </h3>
                <div style={{ display: 'grid', gridTemplateColumns: 'repeat(3, 1fr)', gap: '1rem' }}>
                  {(() => {
                    const allValues = [];
                    for (const m of members) {
                      const data = results.analyzer.getDiagram(m.id, view);
                      allValues.push(...data.map(d => d.val));
                    }
                    const maxVal = Math.max(...allValues.map(v => Math.abs(v)));
                    const maxPos = Math.max(...allValues);
                    const maxNeg = Math.min(...allValues);
                    
                    return (
                      <>
                        <div>
                          <div style={{ fontSize: '0.75rem', color: '#3b82f6', marginBottom: '0.5rem', fontWeight: '600' }}>Maximum Value</div>
                          <div style={{ fontSize: '1.5rem', fontWeight: 'bold', color: '#1e3a8a' }}>{maxVal.toFixed(3)}</div>
                          <div style={{ fontSize: '0.75rem', color: '#6b7280' }}>{view === 'bmd' ? 'kN·m' : 'kN'}</div>
                        </div>
                        <div>
                          <div style={{ fontSize: '0.75rem', color: '#3b82f6', marginBottom: '0.5rem', fontWeight: '600' }}>Max Positive</div>
                          <div style={{ fontSize: '1.5rem', fontWeight: 'bold', color: '#059669' }}>{maxPos.toFixed(3)}</div>
                          <div style={{ fontSize: '0.75rem', color: '#6b7280' }}>{view === 'bmd' ? 'kN·m' : 'kN'}</div>
                        </div>
                        <div>
                          <div style={{ fontSize: '0.75rem', color: '#3b82f6', marginBottom: '0.5rem', fontWeight: '600' }}>Max Negative</div>
                          <div style={{ fontSize: '1.5rem', fontWeight: 'bold', color: '#dc2626' }}>{maxNeg.toFixed(3)}</div>
                          <div style={{ fontSize: '0.75rem', color: '#6b7280' }}>{view === 'bmd' ? 'kN·m' : 'kN'}</div>
                        </div>
                      </>
                    );
                  })()}
                </div>
              </div>

              {/* Individual Member Diagrams */}
              {members.map((member) => {
                const diagramData = results.analyzer.getDiagram(member.id, view, 100);
                const chartData = diagramData.map(d => ({
                  x: d.s.toFixed(2),
                  value: d.val
                }));
                const res = results.memberResults.get(member.id);

                return (
                  <div key={member.id} style={{ ...styles.resultCard, marginBottom: '1.5rem' }}>
                    <h3 style={{ ...styles.sectionTitle, marginBottom: '1rem' }}>
                      Member {member.id + 1} - {view === 'bmd' ? 'Bending Moment' : 'Shear Force'}
                      <span style={{ fontSize: '0.75rem', fontWeight: 'normal', color: '#6b7280', marginLeft: '0.5rem' }}>
                        (Node {member.start + 1} → Node {member.end + 1}, L = {res.L.toFixed(2)}m)
                      </span>
                    </h3>
                    <ResponsiveContainer width="100%" height={250}>
                      <AreaChart data={chartData}>
                        <defs>
                          <linearGradient id={`${view}-grad-${member.id}`} x1="0" y1="0" x2="0" y2="1">
                            <stop offset="0%" stopColor={view === 'bmd' ? '#10b981' : '#3b82f6'} stopOpacity={0.5} />
                            <stop offset="100%" stopColor={view === 'bmd' ? '#10b981' : '#3b82f6'} stopOpacity={0.1} />
                          </linearGradient>
                        </defs>
                        <CartesianGrid strokeDasharray="3 3" stroke="#e5e7eb" />
                        <XAxis 
                          dataKey="x" 
                          stroke="#6b7280" 
                          tick={{ fill: '#6b7280', fontSize: 12 }}
                          label={{ value: 'Distance (m)', position: 'insideBottom', offset: -5, fill: '#6b7280' }}
                        />
                        <YAxis 
                          stroke="#6b7280" 
                          tick={{ fill: '#6b7280', fontSize: 12 }}
                          label={{ value: view === 'bmd' ? 'Moment (kN·m)' : 'Shear (kN)', angle: -90, position: 'insideLeft', fill: '#6b7280' }}
                        />
                        <Tooltip 
                          contentStyle={{ 
                            background: 'white', 
                            border: '1px solid #e5e7eb', 
                            borderRadius: '0.5rem',
                            fontSize: '0.875rem'
                          }}
                          formatter={(value) => [value.toFixed(3), view === 'bmd' ? 'Moment' : 'Shear']}
                          labelFormatter={(label) => `Distance: ${label} m`}
                        />
                        <Area 
                          type="monotone" 
                          dataKey="value" 
                          stroke={view === 'bmd' ? '#10b981' : '#3b82f6'} 
                          strokeWidth={2} 
                          fill={`url(#${view}-grad-${member.id})`} 
                        />
                      </AreaChart>
                    </ResponsiveContainer>
                    <div style={{ marginTop: '1rem', display: 'grid', gridTemplateColumns: 'repeat(4, 1fr)', gap: '0.75rem', fontSize: '0.875rem' }}>
                      <div style={{ padding: '0.5rem', background: '#f9fafb', borderRadius: '0.5rem' }}>
                        <div style={{ fontSize: '0.75rem', color: '#6b7280' }}>Start</div>
                        <div style={{ fontWeight: 'bold', color: '#1e3a8a' }}>
                          {view === 'bmd' ? res.mAB.toFixed(3) : res.vA.toFixed(3)}
                        </div>
                      </div>
                      <div style={{ padding: '0.5rem', background: '#f9fafb', borderRadius: '0.5rem' }}>
                        <div style={{ fontSize: '0.75rem', color: '#6b7280' }}>End</div>
                        <div style={{ fontWeight: 'bold', color: '#1e3a8a' }}>
                          {view === 'bmd' ? res.mBA.toFixed(3) : res.vB.toFixed(3)}
                        </div>
                      </div>
                      <div style={{ padding: '0.5rem', background: '#f9fafb', borderRadius: '0.5rem' }}>
                        <div style={{ fontSize: '0.75rem', color: '#6b7280' }}>Maximum</div>
                        <div style={{ fontWeight: 'bold', color: '#059669' }}>
                          {Math.max(...diagramData.map(d => Math.abs(d.val))).toFixed(3)}
                        </div>
                      </div>
                      <div style={{ padding: '0.5rem', background: '#f9fafb', borderRadius: '0.5rem' }}>
                        <div style={{ fontSize: '0.75rem', color: '#6b7280' }}>Length</div>
                        <div style={{ fontWeight: 'bold', color: '#6b7280' }}>
                          {res.L.toFixed(2)} m
                        </div>
                      </div>
                    </div>
                  </div>
                );
              })}
            </>
          )}

          {/* Results */}
          {results && (
            <>
              <div style={{ display: 'grid', gridTemplateColumns: 'repeat(4, 1fr)', gap: '1rem', marginBottom: '1.5rem' }}>
                <div style={{ ...styles.gradientCard, background: 'linear-gradient(135deg, #dbeafe, #bfdbfe)' }}>
                  <div style={{ fontSize: '0.75rem', color: '#1e40af', marginBottom: '0.5rem' }}>Max Beam Moment</div>
                  <div style={{ fontSize: '1.5rem', fontWeight: 'bold', color: '#1e3a8a' }}>{results.maxMB.toFixed(2)}</div>
                  <div style={{ fontSize: '0.75rem', color: '#3b82f6' }}>kN·m</div>
                </div>
                <div style={{ ...styles.gradientCard, background: 'linear-gradient(135deg, #d1fae5, #a7f3d0)' }}>
                  <div style={{ fontSize: '0.75rem', color: '#065f46', marginBottom: '0.5rem' }}>Max Column Moment</div>
                  <div style={{ fontSize: '1.5rem', fontWeight: 'bold', color: '#064e3b' }}>{results.maxMC.toFixed(2)}</div>
                  <div style={{ fontSize: '0.75rem', color: '#10b981' }}>kN·m</div>
                </div>
                <div style={{ ...styles.gradientCard, background: 'linear-gradient(135deg, #e9d5ff, #d8b4fe)' }}>
                  <div style={{ fontSize: '0.75rem', color: '#6b21a8', marginBottom: '0.5rem' }}>Sway</div>
                  <div style={{ fontSize: '1.5rem', fontWeight: 'bold', color: '#581c87' }}>{(results.swayDelta * 1000).toFixed(2)}</div>
                  <div style={{ fontSize: '0.75rem', color: '#a855f7' }}>mm</div>
                </div>
                <div style={{ ...styles.gradientCard, background: 'linear-gradient(135deg, #fed7aa, #fdba74)' }}>
                  <div style={{ fontSize: '0.75rem', color: '#92400e', marginBottom: '0.5rem' }}>Status</div>
                  <div style={{ fontSize: '1.25rem', fontWeight: 'bold', color: '#78350f' }}>✓ Complete</div>
                </div>
              </div>

              {/* Design Results */}
              <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '1.5rem', marginBottom: '1.5rem' }}>
                {results.beamDes && (
                  <div style={styles.resultCard}>
                    <h3 style={{ ...styles.sectionTitle, marginBottom: '1rem' }}>Beam Design</h3>
                    <div style={{ display: 'flex', flexDirection: 'column', gap: '0.75rem', fontSize: '0.875rem' }}>
                      <div style={{ display: 'flex', justifyContent: 'space-between' }}>
                        <span style={{ color: '#6b7280' }}>Depth:</span>
                        <span style={{ fontWeight: '600' }}>{results.beamDes.D} mm</span>
                      </div>
                      <div style={{ display: 'flex', justifyContent: 'space-between' }}>
                        <span style={{ color: '#6b7280' }}>Steel Area:</span>
                        <span style={{ fontWeight: '600' }}>{results.beamDes.Ast.toFixed(0)} mm²</span>
                      </div>
                      <div style={{ display: 'flex', justifyContent: 'space-between' }}>
                        <span style={{ color: '#6b7280' }}>Type:</span>
                        <span style={{ fontWeight: '600', color: '#3b82f6' }}>{results.beamDes.type}</span>
                      </div>
                    </div>
                  </div>
                )}
                {results.colDes && (
                  <div style={styles.resultCard}>
                    <h3 style={{ ...styles.sectionTitle, marginBottom: '1rem' }}>Column Design</h3>
                    <div style={{ display: 'flex', flexDirection: 'column', gap: '0.75rem', fontSize: '0.875rem' }}>
                      <div style={{ display: 'flex', justifyContent: 'space-between' }}>
                        <span style={{ color: '#6b7280' }}>Size:</span>
                        <span style={{ fontWeight: '600' }}>{results.colDes.size} mm</span>
                      </div>
                      <div style={{ display: 'flex', justifyContent: 'space-between' }}>
                        <span style={{ color: '#6b7280' }}>Steel %:</span>
                        <span style={{ fontWeight: '600' }}>{results.colDes.pct.toFixed(2)}%</span>
                      </div>
                      <div style={{ display: 'flex', justifyContent: 'space-between' }}>
                        <span style={{ color: '#6b7280' }}>Utilization:</span>
                        <span style={{ fontWeight: '600', color: '#10b981' }}>{results.colDes.util.toFixed(1)}%</span>
                      </div>
                    </div>
                  </div>
                )}
              </div>

              {/* Member Results */}
              <div style={styles.resultCard}>
                <h3 style={{ ...styles.sectionTitle, marginBottom: '1rem' }}>Member Forces</h3>
                <div style={{ overflowX: 'auto' }}>
                  <table style={{ width: '100%', fontSize: '0.875rem' }}>
                    <thead>
                      <tr style={{ borderBottom: '2px solid #e5e7eb' }}>
                        <th style={{ padding: '0.75rem', textAlign: 'left', fontWeight: '600', color: '#374151' }}>Member</th>
                        <th style={{ padding: '0.75rem', textAlign: 'right', fontWeight: '600', color: '#374151' }}>M_AB</th>
                        <th style={{ padding: '0.75rem', textAlign: 'right', fontWeight: '600', color: '#374151' }}>M_BA</th>
                        <th style={{ padding: '0.75rem', textAlign: 'right', fontWeight: '600', color: '#374151' }}>V_A</th>
                        <th style={{ padding: '0.75rem', textAlign: 'right', fontWeight: '600', color: '#374151' }}>V_B</th>
                      </tr>
                    </thead>
                    <tbody style={{ fontFamily: 'monospace', fontSize: '0.75rem' }}>
                      {Array.from(results.memberResults.entries()).map(([id, res]) => (
                        <tr key={id} style={{ borderBottom: '1px solid #f3f4f6' }}>
                          <td style={{ padding: '0.75rem', color: '#3b82f6', fontWeight: '600' }}>{id + 1}</td>
                          <td style={{ padding: '0.75rem', textAlign: 'right' }}>{res.mAB.toFixed(3)}</td>
                          <td style={{ padding: '0.75rem', textAlign: 'right' }}>{res.mBA.toFixed(3)}</td>
                          <td style={{ padding: '0.75rem', textAlign: 'right' }}>{res.vA.toFixed(3)}</td>
                          <td style={{ padding: '0.75rem', textAlign: 'right' }}>{res.vB.toFixed(3)}</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              </div>
            </>
          )}
        </main>
      </div>

      <style>{`
        input:focus, select:focus {
          border-color: #3b82f6 !important;
          box-shadow: 0 0 0 3px rgba(59, 130, 246, 0.1);
        }
      `}</style>
    </div>
  );
}