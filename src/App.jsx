import { useState } from 'react';
import BeamAnalyzer from './beam-analyzer';
import FrameAnalyzer from './frame-analyzer';

function App() {
  const [mode, setMode] = useState('beam');
  
  return (
    <div>
      <div style={{ 
        padding: '10px 16px', 
        background: '#1e293b',
        display: 'flex',
        gap: '8px',
        flexWrap: 'wrap'
      }}>
        <button 
          onClick={() => setMode('beam')}
          style={{ 
            padding: '8px 16px',
            background: mode === 'beam' ? '#3b82f6' : 'rgba(255,255,255,0.1)',
            color: 'white',
            border: 'none',
            borderRadius: '6px',
            cursor: 'pointer',
            fontSize: '0.875rem',
            fontWeight: '500',
            flex: '1',
            minWidth: '120px',
          }}
        >
          Beam Analyzer
        </button>
        <button 
          onClick={() => setMode('frame')}
          style={{ 
            padding: '8px 16px',
            background: mode === 'frame' ? '#3b82f6' : 'rgba(255,255,255,0.1)',
            color: 'white',
            border: 'none',
            borderRadius: '6px',
            cursor: 'pointer',
            fontSize: '0.875rem',
            fontWeight: '500',
            flex: '1',
            minWidth: '120px',
          }}
        >
          Frame Analyzer
        </button>
      </div>
      {mode === 'beam' ? <BeamAnalyzer /> : <FrameAnalyzer />}
    </div>
  );
}

export default App;
