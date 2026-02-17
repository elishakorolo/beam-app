import { useState } from 'react';
import BeamAnalyzer from './beam-analyzer';
import FrameAnalyzer from './frame-analyzer';

function App() {
  const [mode, setMode] = useState('beam');
  
  return (
    <div>
      <div style={{ padding: '20px', background: '#f0f0f0' }}>
        <button 
          onClick={() => setMode('beam')}
          style={{ 
            padding: '10px 20px', 
            marginRight: '10px',
            background: mode === 'beam' ? '#007bff' : '#ccc',
            color: 'white',
            border: 'none',
            borderRadius: '5px',
            cursor: 'pointer'
          }}
        >
          Beam Analyzer
        </button>
        <button 
          onClick={() => setMode('frame')}
          style={{ 
            padding: '10px 20px',
            background: mode === 'frame' ? '#007bff' : '#ccc',
            color: 'white',
            border: 'none',
            borderRadius: '5px',
            cursor: 'pointer'
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