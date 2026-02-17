import React, { useState } from 'react';
import { AreaChart, Area, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer } from 'recharts';

// ============================================================================
// STRUCTURAL ANALYSIS ENGINE
// ============================================================================

class StructuralAnalyzer {
  constructor(nodes, E, I_base, loads) {
    this.nodes = nodes.sort((a, b) => a.position - b.position);
    this.E = E * 1e9;
    this.I_base = I_base * 1e-6;
    this.loads = loads;
    this.elements = [];
    
    for (let i = 0; i < this.nodes.length - 1; i++) {
      const I_mult = this.nodes[i].I || 1;
      const length = this.nodes[i + 1].position - this.nodes[i].position;
      
      if (length <= 0) {
        throw new Error(`Invalid span length between nodes ${i} and ${i + 1}. Check node positions.`);
      }
      
      if (I_mult <= 0 || isNaN(I_mult)) {
        throw new Error(`Invalid I multiplier at node ${i}. Must be positive.`);
      }
      
      this.elements.push({
        startNode: this.nodes[i],
        endNode: this.nodes[i + 1],
        length: length,
        stiffness: this.E * this.I_base * I_mult,
        momentStart: 0,
        momentEnd: 0,
        shearStart: 0,
        shearEnd: 0
      });
    }
  }

  calculateFixedEndMoments(element) {
    let femI = 0, femJ = 0;
    const L = element.length;
    const xi = element.startNode.position;

    for (const load of this.loads) {
      if (load.loadType === 'point') {
        if (xi <= load.position && load.position <= xi + L) {
          const a = load.position - xi;
          const b = L - a;
          femI -= (load.magnitude * a * b * b) / (L * L);
          femJ += (load.magnitude * a * a * b) / (L * L);
        }
      } else if (load.loadType === 'udl' && load.endPosition) {
        const start = Math.max(load.position, xi);
        const end = Math.min(load.endPosition, xi + L);
        if (start < end) {
          const w = load.magnitude;
          const c = end - start;
          if (Math.abs(c - L) < 1e-5) {
            femI -= (w * L * L) / 12;
            femJ += (w * L * L) / 12;
          } else {
            const P = w * c;
            const xc = start - xi + c / 2;
            const bc = L - xc;
            femI -= (P * xc * bc * bc) / (L * L);
            femJ += (P * xc * xc * bc) / (L * L);
          }
        }
      } else if (load.loadType === 'vdl' && load.endPosition && load.endMagnitude !== undefined) {
        const start = Math.max(load.position, xi);
        const end = Math.min(load.endPosition, xi + L);
        if (start < end) {
          const totalLen = load.endPosition - load.position;
          const w1 = load.magnitude;
          const w2 = load.endMagnitude;
          
          const ratioStart = totalLen > 1e-6 ? (start - load.position) / totalLen : 0;
          const wStart = w1 + ratioStart * (w2 - w1);
          const ratioEnd = totalLen > 1e-6 ? (end - load.position) / totalLen : 1;
          const wEnd = w1 + ratioEnd * (w2 - w1);
          
          const a = start - xi;
          const c = end - start;
          
          if (Math.abs(a) < 1e-6 && Math.abs(c - L) < 1e-6) {
            femI -= (3 * wStart + 2 * wEnd) * L * L / 60;
            femJ += (2 * wStart + 3 * wEnd) * L * L / 60;
          } else {
            const Ptotal = (wStart + wEnd) / 2 * c;
            const cBar = Math.abs(wStart + wEnd) > 1e-6 
              ? (wStart + 2 * wEnd) * c / (3 * (wStart + wEnd))
              : c / 2;
            const xc = a + cBar;
            const bc = L - xc;
            femI -= (Ptotal * xc * bc * bc) / (L * L);
            femJ += (Ptotal * xc * xc * bc) / (L * L);
          }
        }
      }
    }

    return [femI, femJ];
  }

  getCantileverMoment(node, side) {
    let moment = 0;
    let targetElem = null;

    if (side === 'left') {
      targetElem = this.elements.find(e => 
        e.endNode === node && e.startNode.supportType === 'free'
      ) || null;
    } else {
      targetElem = this.elements.find(e => 
        e.startNode === node && e.endNode.supportType === 'free'
      ) || null;
    }

    if (!targetElem) return 0;

    const refPoint = side === 'left' ? targetElem.endNode.position : targetElem.startNode.position;

    for (const load of this.loads) {
      if (load.loadType === 'point') {
        if (targetElem.startNode.position <= load.position && 
            load.position <= targetElem.endNode.position) {
          const dist = Math.abs(refPoint - load.position);
          moment += load.magnitude * dist;
        }
      } else if (load.loadType === 'udl' && load.endPosition) {
        const s = Math.max(load.position, targetElem.startNode.position);
        const e = Math.min(load.endPosition, targetElem.endNode.position);
        if (s < e) {
          const length = e - s;
          const centroid = s + length / 2;
          const dist = Math.abs(refPoint - centroid);
          moment += (load.magnitude * length) * dist;
        }
      } else if (load.loadType === 'vdl' && load.endPosition && load.endMagnitude !== undefined) {
        const s = Math.max(load.position, targetElem.startNode.position);
        const e = Math.min(load.endPosition, targetElem.endNode.position);
        if (s < e) {
          const totalLen = load.endPosition - load.position;
          const w1 = load.magnitude;
          const w2 = load.endMagnitude;
          
          const ratioStart = totalLen > 1e-6 ? (s - load.position) / totalLen : 0;
          const wStart = w1 + ratioStart * (w2 - w1);
          const ratioEnd = totalLen > 1e-6 ? (e - load.position) / totalLen : 1;
          const wEnd = w1 + ratioEnd * (w2 - w1);
          
          const length = e - s;
          const Ptotal = (wStart + wEnd) / 2 * length;
          const cBar = Math.abs(wStart + wEnd) > 1e-6 
            ? (wStart + 2 * wEnd) * length / (3 * (wStart + wEnd))
            : length / 2;
          const centroid = s + cBar;
          const dist = Math.abs(refPoint - centroid);
          moment += Ptotal * dist;
        }
      }
    }

    return moment;
  }

  solve() {
    const dofNodes = this.nodes.filter(n => 
      n.supportType === 'pinned' || n.supportType === 'roller'
    );
    const n = dofNodes.length;

    if (n === 0) return;

    const K = Array(n).fill(0).map(() => Array(n).fill(0));
    const F = Array(n).fill(0);
    const nodeIndexMap = new Map(dofNodes.map((node, idx) => [node.id, idx]));

    for (const elem of this.elements) {
      if (elem.startNode.supportType === 'free' || elem.endNode.supportType === 'free') 
        continue;

      const k = 2 * elem.stiffness / elem.length;
      const [femI, femJ] = this.calculateFixedEndMoments(elem);

      const iIdx = nodeIndexMap.get(elem.startNode.id);
      const jIdx = nodeIndexMap.get(elem.endNode.id);

      if (iIdx !== undefined) {
        K[iIdx][iIdx] += 2 * k;
        F[iIdx] -= femI;
        if (jIdx !== undefined) K[iIdx][jIdx] += k;
      }

      if (jIdx !== undefined) {
        K[jIdx][jIdx] += 2 * k;
        F[jIdx] -= femJ;
        if (iIdx !== undefined) K[jIdx][iIdx] += k;
      }
    }

    for (let i = 0; i < dofNodes.length; i++) {
      const mLeft = this.getCantileverMoment(dofNodes[i], 'left');
      const mRight = this.getCantileverMoment(dofNodes[i], 'right');
      F[i] += mLeft + mRight;
    }

    const thetas = this.solveLinearSystem(K, F);
    dofNodes.forEach((node, i) => node.rotation = thetas[i]);

    for (const elem of this.elements) {
      const k = 2 * elem.stiffness / elem.length;
      const [femI, femJ] = this.calculateFixedEndMoments(elem);

      if (elem.startNode.supportType === 'free') {
        elem.momentStart = 0;
        elem.momentEnd = -this.getCantileverMoment(elem.endNode, 'left');
      } else if (elem.endNode.supportType === 'free') {
        elem.momentStart = -this.getCantileverMoment(elem.startNode, 'right');
        elem.momentEnd = 0;
      } else {
        elem.momentStart = femI + k * (2 * elem.startNode.rotation + elem.endNode.rotation);
        elem.momentEnd = femJ + k * (elem.startNode.rotation + 2 * elem.endNode.rotation);
      }

      let totalLoad = 0;
      let momentAboutJ = 0;

      for (const load of this.loads) {
        if (load.loadType === 'point') {
          if (elem.startNode.position <= load.position && 
              load.position <= elem.endNode.position) {
            const distFromJ = elem.endNode.position - load.position;
            momentAboutJ += load.magnitude * distFromJ;
            totalLoad += load.magnitude;
          }
        } else if (load.loadType === 'udl' && load.endPosition) {
          const s = Math.max(load.position, elem.startNode.position);
          const e = Math.min(load.endPosition, elem.endNode.position);
          if (s < e) {
            const P = load.magnitude * (e - s);
            const centroid = s + (e - s) / 2;
            const distFromJ = elem.endNode.position - centroid;
            momentAboutJ += P * distFromJ;
            totalLoad += P;
          }
        } else if (load.loadType === 'vdl' && load.endPosition && load.endMagnitude !== undefined) {
          const s = Math.max(load.position, elem.startNode.position);
          const e = Math.min(load.endPosition, elem.endNode.position);
          if (s < e) {
            const totalLen = load.endPosition - load.position;
            const w1 = load.magnitude;
            const w2 = load.endMagnitude;
            
            const ratioStart = totalLen > 1e-6 ? (s - load.position) / totalLen : 0;
            const wStart = w1 + ratioStart * (w2 - w1);
            const ratioEnd = totalLen > 1e-6 ? (e - load.position) / totalLen : 1;
            const wEnd = w1 + ratioEnd * (w2 - w1);
            
            const length = e - s;
            const P = (wStart + wEnd) / 2 * length;
            const cBar = Math.abs(wStart + wEnd) > 1e-6 
              ? (wStart + 2 * wEnd) * length / (3 * (wStart + wEnd))
              : length / 2;
            const centroid = s + cBar;
            const distFromJ = elem.endNode.position - centroid;
            momentAboutJ += P * distFromJ;
            totalLoad += P;
          }
        }
      }

      elem.shearStart = (momentAboutJ - (elem.momentStart + elem.momentEnd)) / elem.length;
      elem.shearEnd = -(totalLoad - elem.shearStart);
    }

    for (const node of this.nodes) {
      if (node.supportType === 'free') continue;
      let reaction = 0;
      for (const elem of this.elements) {
        if (elem.startNode === node) reaction += elem.shearStart;
        if (elem.endNode === node) reaction -= elem.shearEnd;
      }
      node.reaction = reaction;
    }
  }

  solveLinearSystem(A, b) {
    const n = A.length;
    const augmented = A.map((row, i) => [...row, b[i]]);

    for (let i = 0; i < n; i++) {
      let maxRow = i;
      for (let k = i + 1; k < n; k++) {
        if (Math.abs(augmented[k][i]) > Math.abs(augmented[maxRow][i])) {
          maxRow = k;
        }
      }
      [augmented[i], augmented[maxRow]] = [augmented[maxRow], augmented[i]];

      if (Math.abs(augmented[i][i]) < 1e-10) {
        throw new Error('Unstable structure - check support conditions and span lengths');
      }

      for (let k = i + 1; k < n; k++) {
        const factor = augmented[k][i] / augmented[i][i];
        for (let j = i; j < n + 1; j++) {
          augmented[k][j] -= factor * augmented[i][j];
        }
      }
    }

    const x = Array(n).fill(0);
    for (let i = n - 1; i >= 0; i--) {
      x[i] = augmented[i][n];
      for (let j = i + 1; j < n; j++) {
        x[i] -= augmented[i][j] * x[j];
      }
      
      if (Math.abs(augmented[i][i]) < 1e-10) {
        throw new Error('Singular matrix - structure is unstable');
      }
      
      x[i] /= augmented[i][i];
      
      if (isNaN(x[i])) {
        throw new Error('Calculation error - check input values');
      }
    }

    return x;
  }

  calculateDiagrams(numPoints = 500) {
    const span = this.nodes[this.nodes.length - 1].position - this.nodes[0].position;
    const shearData = [];
    const momentData = [];
    let maxMoment = 0;

    for (let i = 0; i <= numPoints; i++) {
      const x = this.nodes[0].position + (span * i) / numPoints;
      
      let shear = 0;
      let moment = 0;

      for (const node of this.nodes) {
        if (node.position <= x && node.supportType !== 'free') {
          shear += node.reaction;
          moment += node.reaction * (x - node.position);
        }
      }

      for (const load of this.loads) {
        if (load.loadType === 'point' && load.position <= x) {
          shear -= load.magnitude;
          moment -= load.magnitude * (x - load.position);
        } else if (load.loadType === 'udl' && load.endPosition && load.position <= x) {
          const endLoad = Math.min(x, load.endPosition);
          if (endLoad > load.position) {
            const length = endLoad - load.position;
            const P = load.magnitude * length;
            const centroid = load.position + length / 2;
            shear -= P;
            moment -= P * (x - centroid);
          }
        } else if (load.loadType === 'vdl' && load.endPosition && 
                   load.endMagnitude !== undefined && load.position <= x) {
          const endLoad = Math.min(x, load.endPosition);
          if (endLoad > load.position) {
            const length = endLoad - load.position;
            const totalLen = load.endPosition - load.position;
            const w1 = load.magnitude;
            const w2 = load.endMagnitude;
            
            const ratioEnd = totalLen > 1e-6 ? length / totalLen : 1;
            const wEnd = w1 + ratioEnd * (w2 - w1);
            
            const P = (w1 + wEnd) / 2 * length;
            const cBar = Math.abs(w1 + wEnd) > 1e-6 
              ? (w1 + 2 * wEnd) * length / (3 * (w1 + wEnd))
              : length / 2;
            const centroid = load.position + cBar;
            shear -= P;
            moment -= P * (x - centroid);
          }
        }
      }

      shearData.push({ x, value: shear });
      momentData.push({ x, value: moment });
      maxMoment = Math.max(maxMoment, Math.abs(moment));
    }

    return { shear: shearData, moment: momentData, maxMoment };
  }
}

// ============================================================================
// DESIGN CALCULATORS
// ============================================================================

function designConcreteBeam(maxMoment, params) {
  const fck = params.fck || 25;
  const fy = params.fy || 415;
  const b = (params.beamWidth || 0.3) * 1000;
  const cover = params.cover || 40;
  const dia = params.barDia || 16;

  const fcd = 0.67 * fck / 1.5;
  const fyd = fy / 1.15;

  const xuMaxRatio = 0.0035 / (0.0055 + 0.87 * fy / 200000);
  const MuLimRatio = 0.36 * xuMaxRatio * (1 - 0.42 * xuMaxRatio) * fcd;

  const Mdesign = Math.abs(maxMoment) * 1e6;
  const dReq = Math.sqrt(Mdesign / (MuLimRatio * b));
  let D = dReq + cover + dia / 2;
  D = Math.ceil(D / 25) * 25;
  const d = D - cover - dia / 2;

  const MuLim = MuLimRatio * b * d * d;
  let Ast, Asc, sectionType;

  if (Mdesign <= MuLim) {
    const p = Mdesign / (b * d * d * fcd);
    Ast = (b * d * fcd / fyd) * (1 - Math.sqrt(1 - 4.6 * p));
    Asc = 0;
    sectionType = "Singly Reinforced";
  } else {
    const Ast1 = MuLimRatio * b * d * fcd / fyd;
    const Mu2 = Mdesign - MuLim;
    const dPrime = cover + dia;
    Asc = Mu2 / (fyd * (d - dPrime));
    Ast = Ast1 + Asc;
    sectionType = "Doubly Reinforced";
  }

  const AstMin = 0.85 * b * d / fy;
  Ast = Math.max(Ast, AstMin);

  const numBars = Math.ceil(Ast / (Math.PI * dia * dia / 4));
  const providedArea = numBars * (Math.PI * dia * dia / 4);

  return {
    depth: D,
    effDepth: d,
    steelArea: Ast,
    compSteel: Asc,
    sectionType,
    numBars,
    providedArea
  };
}

function designSteelBeam(maxMoment, params) {
  const fy = params.fySteelMpa || 250;
  const fd = fy / 1.1;
  const Mdesign = Math.abs(maxMoment);
  const ZreqCm3 = (Mdesign * 1e6) / fd / 1000;

  const sections = {
    'ISMB 100': 87.4, 'ISMB 125': 139, 'ISMB 150': 197, 'ISMB 175': 267,
    'ISMB 200': 353, 'ISMB 225': 453, 'ISMB 250': 571, 'ISMB 300': 788,
    'ISMB 350': 1051, 'ISMB 400': 1377, 'ISMB 450': 1753, 'ISMB 500': 2194,
    'ISMB 550': 2691, 'ISMB 600': 3277
  };

  let selectedSection = 'Custom (> ISMB 600)';
  let utilization = 100;

  for (const [section, Z] of Object.entries(sections)) {
    if (Z >= ZreqCm3) {
      selectedSection = section;
      utilization = (ZreqCm3 / Z) * 100;
      break;
    }
  }

  return { requiredZ: ZreqCm3, section: selectedSection, utilization };
}

// ============================================================================
// MINIMAL CLEAN STYLES
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
    padding: '0.5rem 1rem',
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
  deleteBtn: {
    color: '#ef4444',
    background: 'transparent',
    border: 'none',
    cursor: 'pointer',
    fontSize: '0.75rem'
  },
  mainContent: {
    flex: 1,
    padding: '2rem',
    overflowY: 'auto'
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
// REACT COMPONENT
// ============================================================================

export default function App() {
  const [totalLength, setTotalLength] = useState(10);
  const [youngModulus, setYoungModulus] = useState(200);
  const [momentInertia, setMomentInertia] = useState(1);
  const [nodes, setNodes] = useState([
    { id: 0, position: 0, supportType: 'pinned', rotation: 0, reaction: 0, I: 1 },
    { id: 1, position: 5, supportType: 'roller', rotation: 0, reaction: 0, I: 1 },
    { id: 2, position: 10, supportType: 'pinned', rotation: 0, reaction: 0, I: 1 }
  ]);
  const [loads, setLoads] = useState([
    { id: 0, loadType: 'point', magnitude: 10, position: 2.5 }
  ]);
  const [designParams, setDesignParams] = useState({
    material: 'concrete',
    fck: 25,
    fy: 415,
    beamWidth: 0.3,
    cover: 40,
    barDia: 16,
    fySteelMpa: 250
  });

  const [results, setResults] = useState(null);
  const [activeTab, setActiveTab] = useState('nodes');

  const addNode = () => {
    const newId = nodes.length;
    const newPos = nodes.length > 0 ? nodes[nodes.length - 1].position + 2 : 0;
    setNodes([...nodes, {
      id: newId,
      position: newPos,
      supportType: 'roller',
      rotation: 0,
      reaction: 0,
      I: 1
    }]);
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

  const analyze = () => {
    try {
      const analyzer = new StructuralAnalyzer([...nodes], youngModulus, momentInertia, loads);
      analyzer.solve();
      const diagrams = analyzer.calculateDiagrams(500);
      
      let designResult;
      if (designParams.material === 'concrete') {
        designResult = designConcreteBeam(diagrams.maxMoment, designParams);
      } else {
        designResult = designSteelBeam(diagrams.maxMoment, designParams);
      }

      setResults({
        nodes: [...nodes],
        shearData: diagrams.shear,
        momentData: diagrams.moment,
        maxMoment: diagrams.maxMoment,
        design: designResult
      });
    } catch (error) {
      console.error('Analysis failed:', error);
      alert('Analysis failed. Please check your inputs.');
    }
  };

  return (
    <div style={styles.container}>
      {/* Header */}
      <header style={styles.header}>
        <div style={{ padding: '1rem 1.5rem', display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
          <div>
            <h1 style={styles.title}>Beam Analyzer</h1>
            <p style={styles.subtitle}>Slope-Deflection Analysis & Design</p>
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
        {/* Sidebar */}
        <aside style={styles.sidebar}>
          <div style={{ display: 'flex', gap: '0.25rem', background: '#f3f4f6', padding: '0.25rem', borderRadius: '0.5rem', marginBottom: '1.5rem' }}>
            {['nodes', 'loads', 'design'].map((tab) => (
              <button
                key={tab}
                onClick={() => setActiveTab(tab)}
                style={{
                  ...styles.tab,
                  flex: 1,
                  ...(activeTab === tab ? styles.tabActive : {})
                }}
              >
                {tab.charAt(0).toUpperCase() + tab.slice(1)}
              </button>
            ))}
          </div>

          {/* Global Parameters */}
          <section style={{ marginBottom: '1.5rem' }}>
            <h2 style={styles.sectionTitle}>Geometry</h2>
            <div style={{ display: 'flex', flexDirection: 'column', gap: '1rem' }}>
              <div>
                <label style={styles.label}>E (GPa)</label>
                <input
                  type="number"
                  value={youngModulus}
                  onChange={(e) => setYoungModulus(parseFloat(e.target.value) || 200)}
                  style={styles.input}
                />
              </div>
              <div>
                <label style={styles.label}>I (×10⁻⁶ m⁴)</label>
                <input
                  type="number"
                  value={momentInertia}
                  onChange={(e) => setMomentInertia(parseFloat(e.target.value) || 1)}
                  style={styles.input}
                />
              </div>
            </div>
          </section>

          {/* Nodes Tab */}
          {activeTab === 'nodes' && (
            <section>
              <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1rem' }}>
                <h2 style={styles.sectionTitle}>Supports</h2>
                <button onClick={addNode} style={styles.button}>+ Add</button>
              </div>
              <div style={{ display: 'flex', flexDirection: 'column', gap: '0.75rem', maxHeight: '400px', overflowY: 'auto' }}>
                {nodes.map((node, idx) => (
                  <div key={node.id} style={styles.card}>
                    <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '0.5rem' }}>
                      <span style={{ fontSize: '0.75rem', color: '#3b82f6', fontWeight: '600' }}>Node {idx}</span>
                      {nodes.length > 2 && (
                        <button onClick={() => deleteNode(idx)} style={styles.deleteBtn}>✕</button>
                      )}
                    </div>
                    <div style={{ display: 'flex', flexDirection: 'column', gap: '0.5rem' }}>
                      <div>
                        <label style={styles.label}>Position (m)</label>
                        <input
                          type="number"
                          value={node.position}
                          onChange={(e) => updateNode(idx, 'position', parseFloat(e.target.value) || 0)}
                          style={styles.input}
                        />
                      </div>
                      <div>
                        <label style={styles.label}>Support Type</label>
                        <select
                          value={node.supportType}
                          onChange={(e) => updateNode(idx, 'supportType', e.target.value)}
                          style={styles.select}
                        >
                          <option value="free">Free</option>
                          <option value="pinned">Pinned</option>
                          <option value="roller">Roller</option>
                          <option value="fixed">Fixed</option>
                        </select>
                      </div>
                      {idx < nodes.length - 1 && (
                        <div>
                          <label style={styles.label}>Span I Multiplier</label>
                          <input
                            type="number"
                            value={node.I}
                            onChange={(e) => updateNode(idx, 'I', parseFloat(e.target.value) || 1)}
                            style={styles.input}
                            step="0.1"
                          />
                        </div>
                      )}
                    </div>
                  </div>
                ))}
              </div>
            </section>
          )}

          {/* Loads Tab */}
          {activeTab === 'loads' && (
            <section>
              <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1rem' }}>
                <h2 style={styles.sectionTitle}>Loads</h2>
                <button
                  onClick={() => setLoads([...loads, {
                    id: loads.length,
                    loadType: 'point',
                    magnitude: 10,
                    position: 5
                  }])}
                  style={styles.button}
                >
                  + Add
                </button>
              </div>
              <div style={{ display: 'flex', flexDirection: 'column', gap: '0.75rem', maxHeight: '400px', overflowY: 'auto' }}>
                {loads.map((load, idx) => (
                  <div key={load.id} style={styles.card}>
                    <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '0.5rem' }}>
                      <span style={{ fontSize: '0.75rem', color: '#3b82f6', fontWeight: '600' }}>Load {idx + 1}</span>
                      <button
                        onClick={() => setLoads(loads.filter((_, i) => i !== idx))}
                        style={styles.deleteBtn}
                      >
                        ✕
                      </button>
                    </div>
                    <div style={{ display: 'flex', flexDirection: 'column', gap: '0.5rem' }}>
                      <select
                        value={load.loadType}
                        onChange={(e) => {
                          const newLoads = [...loads];
                          newLoads[idx].loadType = e.target.value;
                          if (e.target.value === 'udl' || e.target.value === 'vdl') {
                            newLoads[idx].endPosition = (load.position || 0) + 2;
                          }
                          if (e.target.value === 'vdl') {
                            newLoads[idx].endMagnitude = 20;
                          }
                          setLoads(newLoads);
                        }}
                        style={styles.select}
                      >
                        <option value="point">Point Load</option>
                        <option value="udl">UDL</option>
                        <option value="vdl">VDL</option>
                      </select>
                      <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '0.5rem' }}>
                        <div>
                          <label style={styles.label}>Mag (kN{load.loadType !== 'point' ? '/m' : ''})</label>
                          <input
                            type="number"
                            value={load.magnitude}
                            onChange={(e) => {
                              const newLoads = [...loads];
                              newLoads[idx].magnitude = parseFloat(e.target.value) || 0;
                              setLoads(newLoads);
                            }}
                            style={styles.input}
                          />
                        </div>
                        <div>
                          <label style={styles.label}>{load.loadType === 'point' ? 'Pos' : 'Start'} (m)</label>
                          <input
                            type="number"
                            value={load.position}
                            onChange={(e) => {
                              const newLoads = [...loads];
                              newLoads[idx].position = parseFloat(e.target.value) || 0;
                              setLoads(newLoads);
                            }}
                            style={styles.input}
                          />
                        </div>
                      </div>
                      {(load.loadType === 'udl' || load.loadType === 'vdl') && (
                        <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '0.5rem' }}>
                          <div>
                            <label style={styles.label}>End (m)</label>
                            <input
                              type="number"
                              value={load.endPosition || 0}
                              onChange={(e) => {
                                const newLoads = [...loads];
                                newLoads[idx].endPosition = parseFloat(e.target.value) || 0;
                                setLoads(newLoads);
                              }}
                              style={styles.input}
                            />
                          </div>
                          {load.loadType === 'vdl' && (
                            <div>
                              <label style={styles.label}>End Mag (kN/m)</label>
                              <input
                                type="number"
                                value={load.endMagnitude || 0}
                                onChange={(e) => {
                                  const newLoads = [...loads];
                                  newLoads[idx].endMagnitude = parseFloat(e.target.value) || 0;
                                  setLoads(newLoads);
                                }}
                                style={styles.input}
                              />
                            </div>
                          )}
                        </div>
                      )}
                    </div>
                  </div>
                ))}
              </div>
            </section>
          )}

          {/* Design Tab */}
          {activeTab === 'design' && (
            <section>
              <h2 style={styles.sectionTitle}>Design Parameters</h2>
              <div style={{ display: 'flex', flexDirection: 'column', gap: '1rem' }}>
                <div>
                  <label style={styles.label}>Material</label>
                  <select
                    value={designParams.material}
                    onChange={(e) => setDesignParams({ ...designParams, material: e.target.value })}
                    style={styles.select}
                  >
                    <option value="concrete">Concrete</option>
                    <option value="steel">Steel</option>
                  </select>
                </div>
                {designParams.material === 'concrete' ? (
                  <>
                    <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '0.75rem' }}>
                      <div>
                        <label style={styles.label}>fck (MPa)</label>
                        <input
                          type="number"
                          value={designParams.fck}
                          onChange={(e) => setDesignParams({ ...designParams, fck: parseFloat(e.target.value) || 25 })}
                          style={styles.input}
                        />
                      </div>
                      <div>
                        <label style={styles.label}>fy (MPa)</label>
                        <input
                          type="number"
                          value={designParams.fy}
                          onChange={(e) => setDesignParams({ ...designParams, fy: parseFloat(e.target.value) || 415 })}
                          style={styles.input}
                        />
                      </div>
                      <div>
                        <label style={styles.label}>Width (m)</label>
                        <input
                          type="number"
                          value={designParams.beamWidth}
                          onChange={(e) => setDesignParams({ ...designParams, beamWidth: parseFloat(e.target.value) || 0.3 })}
                          style={styles.input}
                        />
                      </div>
                      <div>
                        <label style={styles.label}>Cover (mm)</label>
                        <input
                          type="number"
                          value={designParams.cover}
                          onChange={(e) => setDesignParams({ ...designParams, cover: parseFloat(e.target.value) || 40 })}
                          style={styles.input}
                        />
                      </div>
                    </div>
                  </>
                ) : (
                  <div>
                    <label style={styles.label}>fy (MPa)</label>
                    <input
                      type="number"
                      value={designParams.fySteelMpa}
                      onChange={(e) => setDesignParams({ ...designParams, fySteelMpa: parseFloat(e.target.value) || 250 })}
                      style={styles.input}
                    />
                  </div>
                )}
              </div>
            </section>
          )}
        </aside>

        {/* Main Content */}
        <main style={styles.mainContent}>
          {results ? (
            <div style={{ display: 'flex', flexDirection: 'column', gap: '2rem' }}>
              {/* Summary Cards */}
              <div style={{ display: 'grid', gridTemplateColumns: 'repeat(3, 1fr)', gap: '1rem' }}>
                <div style={{ ...styles.gradientCard, background: 'linear-gradient(135deg, #dbeafe, #bfdbfe)' }}>
                  <div style={{ fontSize: '0.75rem', color: '#1e40af', marginBottom: '0.5rem' }}>Max Moment</div>
                  <div style={{ fontSize: '1.5rem', fontWeight: 'bold', color: '#1e3a8a' }}>{results.maxMoment.toFixed(3)}</div>
                  <div style={{ fontSize: '0.75rem', color: '#3b82f6' }}>kN·m</div>
                </div>
                <div style={{ ...styles.gradientCard, background: 'linear-gradient(135deg, #d1fae5, #a7f3d0)' }}>
                  <div style={{ fontSize: '0.75rem', color: '#065f46', marginBottom: '0.5rem' }}>Total Reaction</div>
                  <div style={{ fontSize: '1.5rem', fontWeight: 'bold', color: '#064e3b' }}>
                    {results.nodes.reduce((sum, n) => sum + Math.abs(n.reaction), 0).toFixed(3)}
                  </div>
                  <div style={{ fontSize: '0.75rem', color: '#10b981' }}>kN</div>
                </div>
                <div style={{ ...styles.gradientCard, background: 'linear-gradient(135deg, #fed7aa, #fdba74)' }}>
                  <div style={{ fontSize: '0.75rem', color: '#92400e', marginBottom: '0.5rem' }}>Status</div>
                  <div style={{ fontSize: '1.25rem', fontWeight: 'bold', color: '#78350f' }}>✓ Complete</div>
                </div>
              </div>

              {/* Design Results */}
              <div style={styles.resultCard}>
                <h3 style={{ ...styles.sectionTitle, marginBottom: '1rem' }}>Design Output</h3>
                {designParams.material === 'concrete' ? (
                  <div style={{ display: 'grid', gridTemplateColumns: 'repeat(4, 1fr)', gap: '1rem' }}>
                    <div>
                      <div style={{ fontSize: '0.75rem', color: '#6b7280' }}>Depth</div>
                      <div style={{ fontSize: '1.5rem', fontWeight: 'bold', color: '#1e3a8a' }}>{results.design.depth}</div>
                      <div style={{ fontSize: '0.75rem', color: '#6b7280' }}>mm</div>
                    </div>
                    <div>
                      <div style={{ fontSize: '0.75rem', color: '#6b7280' }}>Steel Area</div>
                      <div style={{ fontSize: '1.5rem', fontWeight: 'bold', color: '#059669' }}>{results.design.steelArea.toFixed(0)}</div>
                      <div style={{ fontSize: '0.75rem', color: '#6b7280' }}>mm²</div>
                    </div>
                    <div>
                      <div style={{ fontSize: '0.75rem', color: '#6b7280' }}>Bars</div>
                      <div style={{ fontSize: '1.5rem', fontWeight: 'bold', color: '#3b82f6' }}>{results.design.numBars}×{designParams.barDia}Ø</div>
                    </div>
                    <div>
                      <div style={{ fontSize: '0.75rem', color: '#6b7280' }}>Type</div>
                      <div style={{ fontSize: '0.875rem', fontWeight: 'bold', color: '#6b7280', marginTop: '0.5rem' }}>{results.design.sectionType}</div>
                    </div>
                  </div>
                ) : (
                  <div style={{ display: 'grid', gridTemplateColumns: 'repeat(3, 1fr)', gap: '1rem' }}>
                    <div>
                      <div style={{ fontSize: '0.75rem', color: '#6b7280' }}>Required Z</div>
                      <div style={{ fontSize: '1.5rem', fontWeight: 'bold', color: '#1e3a8a' }}>{results.design.requiredZ.toFixed(1)}</div>
                      <div style={{ fontSize: '0.75rem', color: '#6b7280' }}>cm³</div>
                    </div>
                    <div>
                      <div style={{ fontSize: '0.75rem', color: '#6b7280' }}>Section</div>
                      <div style={{ fontSize: '1.25rem', fontWeight: 'bold', color: '#3b82f6', marginTop: '0.5rem' }}>{results.design.section}</div>
                    </div>
                    <div>
                      <div style={{ fontSize: '0.75rem', color: '#6b7280' }}>Utilization</div>
                      <div style={{ fontSize: '1.5rem', fontWeight: 'bold', color: '#059669' }}>{results.design.utilization.toFixed(1)}%</div>
                    </div>
                  </div>
                )}
              </div>

              {/* Diagrams */}
              <div style={styles.resultCard}>
                <h3 style={{ ...styles.sectionTitle, marginBottom: '1rem' }}>Shear Force Diagram</h3>
                <ResponsiveContainer width="100%" height={250}>
                  <AreaChart data={results.shearData}>
                    <defs>
                      <linearGradient id="shearGrad" x1="0" y1="0" x2="0" y2="1">
                        <stop offset="0%" stopColor="#3b82f6" stopOpacity={0.5} />
                        <stop offset="100%" stopColor="#3b82f6" stopOpacity={0.1} />
                      </linearGradient>
                    </defs>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e5e7eb" />
                    <XAxis dataKey="x" stroke="#6b7280" tick={{ fill: '#6b7280', fontSize: 12 }} />
                    <YAxis stroke="#6b7280" tick={{ fill: '#6b7280', fontSize: 12 }} />
                    <Tooltip contentStyle={{ background: 'white', border: '1px solid #e5e7eb', borderRadius: '0.5rem' }} />
                    <Area type="monotone" dataKey="value" stroke="#3b82f6" strokeWidth={2} fill="url(#shearGrad)" />
                  </AreaChart>
                </ResponsiveContainer>
              </div>

              <div style={styles.resultCard}>
                <h3 style={{ ...styles.sectionTitle, marginBottom: '1rem' }}>Bending Moment Diagram</h3>
                <ResponsiveContainer width="100%" height={250}>
                  <AreaChart data={results.momentData}>
                    <defs>
                      <linearGradient id="momentGrad" x1="0" y1="0" x2="0" y2="1">
                        <stop offset="0%" stopColor="#10b981" stopOpacity={0.5} />
                        <stop offset="100%" stopColor="#10b981" stopOpacity={0.1} />
                      </linearGradient>
                    </defs>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e5e7eb" />
                    <XAxis dataKey="x" stroke="#6b7280" tick={{ fill: '#6b7280', fontSize: 12 }} />
                    <YAxis stroke="#6b7280" tick={{ fill: '#6b7280', fontSize: 12 }} />
                    <Tooltip contentStyle={{ background: 'white', border: '1px solid #e5e7eb', borderRadius: '0.5rem' }} />
                    <Area type="monotone" dataKey="value" stroke="#10b981" strokeWidth={2} fill="url(#momentGrad)" />
                  </AreaChart>
                </ResponsiveContainer>
              </div>

              {/* Results Tables */}
              <div style={styles.resultCard}>
                <h3 style={{ ...styles.sectionTitle, marginBottom: '1rem' }}>Support End Moments</h3>
                <div style={{ overflowX: 'auto' }}>
                  <table style={{ width: '100%', fontSize: '0.875rem', borderCollapse: 'collapse' }}>
                    <thead>
                      <tr style={{ borderBottom: '2px solid #e5e7eb' }}>
                        <th style={{ padding: '0.75rem', textAlign: 'left', fontWeight: '600', color: '#374151' }}>Node</th>
                        <th style={{ padding: '0.75rem', textAlign: 'right', fontWeight: '600', color: '#374151' }}>Position</th>
                        <th style={{ padding: '0.75rem', textAlign: 'right', fontWeight: '600', color: '#374151' }}>Type</th>
                        <th style={{ padding: '0.75rem', textAlign: 'right', fontWeight: '600', color: '#374151' }}>M_left</th>
                        <th style={{ padding: '0.75rem', textAlign: 'right', fontWeight: '600', color: '#374151' }}>M_right</th>
                      </tr>
                    </thead>
                    <tbody style={{ fontFamily: 'monospace', fontSize: '0.75rem' }}>
                      {(() => {
                        const analyzer = new StructuralAnalyzer([...nodes], youngModulus, momentInertia, loads);
                        analyzer.solve();
                        
                        return analyzer.nodes.map((node, idx) => {
                          const leftSpan = analyzer.elements.find(e => e.endNode.id === node.id);
                          const momentFromLeft = leftSpan ? leftSpan.momentEnd : null;
                          const rightSpan = analyzer.elements.find(e => e.startNode.id === node.id);
                          const momentFromRight = rightSpan ? rightSpan.momentStart : null;
                          
                          return (
                            <tr key={node.id} style={{ borderBottom: '1px solid #f3f4f6' }}>
                              <td style={{ padding: '0.75rem', color: '#3b82f6', fontWeight: '600' }}>{idx}</td>
                              <td style={{ padding: '0.75rem', textAlign: 'right' }}>{node.position.toFixed(3)}</td>
                              <td style={{ padding: '0.75rem', textAlign: 'right', textTransform: 'capitalize' }}>{node.supportType}</td>
                              <td style={{ padding: '0.75rem', textAlign: 'right' }}>
                                {momentFromLeft !== null ? momentFromLeft.toFixed(3) : '-'}
                              </td>
                              <td style={{ padding: '0.75rem', textAlign: 'right' }}>
                                {momentFromRight !== null ? momentFromRight.toFixed(3) : '-'}
                              </td>
                            </tr>
                          );
                        });
                      })()}
                    </tbody>
                  </table>
                </div>
              </div>

              <div style={styles.resultCard}>
                <h3 style={{ ...styles.sectionTitle, marginBottom: '1rem' }}>Span Details</h3>
                <div style={{ overflowX: 'auto' }}>
                  <table style={{ width: '100%', fontSize: '0.875rem', borderCollapse: 'collapse' }}>
                    <thead>
                      <tr style={{ borderBottom: '2px solid #e5e7eb' }}>
                        <th style={{ padding: '0.75rem', textAlign: 'left', fontWeight: '600', color: '#374151' }}>Span</th>
                        <th style={{ padding: '0.75rem', textAlign: 'center', fontWeight: '600', color: '#374151' }}>Nodes</th>
                        <th style={{ padding: '0.75rem', textAlign: 'right', fontWeight: '600', color: '#374151' }}>Length</th>
                        <th style={{ padding: '0.75rem', textAlign: 'right', fontWeight: '600', color: '#374151' }}>M_start</th>
                        <th style={{ padding: '0.75rem', textAlign: 'right', fontWeight: '600', color: '#374151' }}>M_end</th>
                        <th style={{ padding: '0.75rem', textAlign: 'right', fontWeight: '600', color: '#374151' }}>V_start</th>
                        <th style={{ padding: '0.75rem', textAlign: 'right', fontWeight: '600', color: '#374151' }}>V_end</th>
                      </tr>
                    </thead>
                    <tbody style={{ fontFamily: 'monospace', fontSize: '0.75rem' }}>
                      {(() => {
                        const analyzer = new StructuralAnalyzer([...nodes], youngModulus, momentInertia, loads);
                        analyzer.solve();
                        
                        return analyzer.elements.map((elem, idx) => (
                          <tr key={idx} style={{ borderBottom: '1px solid #f3f4f6' }}>
                            <td style={{ padding: '0.75rem', color: '#3b82f6', fontWeight: '600' }}>{idx + 1}</td>
                            <td style={{ padding: '0.75rem', textAlign: 'center' }}>
                              {elem.startNode.id} → {elem.endNode.id}
                            </td>
                            <td style={{ padding: '0.75rem', textAlign: 'right' }}>{elem.length.toFixed(3)}</td>
                            <td style={{ padding: '0.75rem', textAlign: 'right' }}>{elem.momentStart.toFixed(3)}</td>
                            <td style={{ padding: '0.75rem', textAlign: 'right' }}>{elem.momentEnd.toFixed(3)}</td>
                            <td style={{ padding: '0.75rem', textAlign: 'right' }}>{elem.shearStart.toFixed(3)}</td>
                            <td style={{ padding: '0.75rem', textAlign: 'right' }}>{elem.shearEnd.toFixed(3)}</td>
                          </tr>
                        ));
                      })()}
                    </tbody>
                  </table>
                </div>
              </div>

              <div style={styles.resultCard}>
                <h3 style={{ ...styles.sectionTitle, marginBottom: '1rem' }}>Support Reactions</h3>
                <div style={{ display: 'flex', flexDirection: 'column', gap: '0.5rem' }}>
                  {results.nodes.filter(n => n.supportType !== 'free').map((node) => (
                    <div key={node.id} style={{ display: 'flex', justifyContent: 'space-between', padding: '0.75rem', background: '#f9fafb', borderRadius: '0.5rem', border: '1px solid #e5e7eb' }}>
                      <span style={{ fontSize: '0.875rem', color: '#6b7280' }}>Node {node.id} ({node.position.toFixed(2)}m)</span>
                      <span style={{ fontSize: '1.125rem', fontWeight: 'bold', color: '#3b82f6' }}>{node.reaction.toFixed(3)} kN</span>
                    </div>
                  ))}
                </div>
              </div>
            </div>
          ) : (
            <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'center', height: '400px' }}>
              <p style={{ color: '#6b7280', fontSize: '1.125rem' }}>
                Configure parameters and click Analyze
              </p>
            </div>
          )}
        </main>
      </div>

      <style>{`
        input:focus, select:focus {
          border-color: #3b82f6 !important;
          box-shadow: 0 0 0 3px rgba(59, 130, 246, 0.1);
        }
        button:hover {
          opacity: 0.9;
        }
      `}</style>
    </div>
  );
}