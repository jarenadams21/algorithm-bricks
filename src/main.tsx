import { StrictMode } from 'react';
import { createRoot } from 'react-dom/client';
import './index.css';
import Dijkstra from './Dijkstra.tsx';

createRoot(document.getElementById('root')!).render(
  <StrictMode>
    <Dijkstra />
  </StrictMode>,
)
