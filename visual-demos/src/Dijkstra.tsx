import React, { useState, useRef, MouseEvent, ChangeEvent, JSX } from "react";

// ------------------------- Types -------------------------
interface NodeType {
  id: number;
  x: number;
  y: number;
}

interface EdgeType {
  id: number;      // ID for the edge itself
  from: number;    // node ID
  to: number;      // node ID
  weight: number;
}

/** adjacencyList[nodeID] = array of { node, weight } for each neighbor */
type AdjacencyList = Record<number, { node: number; weight: number }[]>;
type Distances = Record<number, number>;

// ------------------------- Logic -------------------------

/**
 * Build an undirected adjacency list so that edges go both ways,
 * from -> to and to -> from.
 *
 * If you prefer directed edges, remove the second push() line.
 */
function buildAdjacencyList(nodes: NodeType[], edges: EdgeType[]): AdjacencyList {
  const adjList: AdjacencyList = {};
  // Initialize
  nodes.forEach((n) => {
    adjList[n.id] = [];
  });

  // Populate
  edges.forEach((edge) => {
    adjList[edge.from].push({ node: edge.to, weight: edge.weight });
    // For undirected, also add the reverse direction:
    adjList[edge.to].push({ node: edge.from, weight: edge.weight });
  });

  return adjList;
}

/**
 * Standard Dijkstra using a simple array as a priority queue (O(VE) worst-case).
 * If you need better performance on large graphs, you'd use a real PQ/heap structure.
 */
function dijkstra(adjList: AdjacencyList, startId: number): Distances {
  const distances: Distances = {};
  const visited: Record<number, boolean> = {};
  // Priority queue: array of (node, dist)
  const pq: { node: number; dist: number }[] = [];

  // Initialize
  for (const nodeIdStr in adjList) {
    const nodeId = parseInt(nodeIdStr, 10);
    distances[nodeId] = Infinity;
    visited[nodeId] = false;
  }
  // distance of the start node is 0
  distances[startId] = 0;
  pq.push({ node: startId, dist: 0 });

  // Process
  while (pq.length > 0) {
    // Sort ascending by distance
    pq.sort((a, b) => a.dist - b.dist);
    const top = pq.shift();
    if (!top) break;

    const { node: current, dist: currentDist } = top;
    if (visited[current]) continue;
    visited[current] = true;

    // Update neighbors
    for (const neighbor of adjList[current]) {
      const newDist = currentDist + neighbor.weight;
      if (newDist < distances[neighbor.node]) {
        distances[neighbor.node] = newDist;
        pq.push({ node: neighbor.node, dist: newDist });
      }
    }
  }

  return distances;
}

// ------------------ Demo Component ------------------

/**
 * Default nodes: A small connected example so you can see it working.
 * We'll label them as 1,2,3,4 for clarity.
 */
const DEFAULT_NODES: NodeType[] = [
  { id: 1, x: 100, y: 100 },
  { id: 2, x: 300, y: 100 },
  { id: 3, x: 100, y: 300 },
  { id: 4, x: 300, y: 300 },
];

/**
 * We'll give them edges so everything is connected in an undirected manner.
 * If you remove the "reverse" adjacency, you'd see that for a directed version,
 * some nodes might have Infinity if there's no path from the start to them.
 */
const DEFAULT_EDGES: EdgeType[] = [
  { id: 101, from: 1, to: 2, weight: 5 },
  { id: 102, from: 2, to: 4, weight: 2 },
  { id: 103, from: 1, to: 3, weight: 6 },
  { id: 104, from: 3, to: 4, weight: 10 },
];

// The UI has three modes
type Mode = "addNode" | "addEdge" | "selectStart";

export default function DijkstraDemo(): JSX.Element {
  // Our working graph
  const [nodes, setNodes] = useState<NodeType[]>(DEFAULT_NODES);
  const [edges, setEdges] = useState<EdgeType[]>(DEFAULT_EDGES);

  // We keep a small integer ID for each *newly created node*:
  const nextNodeIdRef = useRef(5); // after 1,2,3,4

  // UI mode
  const [mode, setMode] = useState<Mode>("addNode");
  // For "addEdge" mode, store the first node we clicked
  const [pendingFrom, setPendingFrom] = useState<number | null>(null);

  // Edge weight to use for new edges
  const [edgeWeight, setEdgeWeight] = useState<number>(1);

  // Dijkstra
  const [startNode, setStartNode] = useState<number | null>(null);
  const [distances, setDistances] = useState<Distances>({});

  /**
   * Add Node: click on empty space in the <svg>.
   */
  const handleSvgClick = (e: MouseEvent<SVGSVGElement>) => {
    if (mode !== "addNode") return;
    const { offsetX, offsetY } = e.nativeEvent;

    const newNode: NodeType = {
      id: nextNodeIdRef.current,
      x: offsetX,
      y: offsetY,
    };
    nextNodeIdRef.current++;
    setNodes((prev) => [...prev, newNode]);
  };

  /**
   * Click on a node (circle).
   * Behavior depends on mode:
   *  - addEdge: set "from" node or create edge
   *  - selectStart: set Dijkstra start node
   */
  const handleNodeClick = (nodeId: number, e: MouseEvent<SVGCircleElement>) => {
    e.stopPropagation(); // so we don't also get the SVG's onClick

    if (mode === "addEdge") {
      if (pendingFrom == null) {
        setPendingFrom(nodeId);
      } else {
        // create an edge from pendingFrom -> nodeId
        if (pendingFrom !== nodeId) {
          const newEdge: EdgeType = {
            id: Date.now(), // or your own ID logic
            from: pendingFrom,
            to: nodeId,
            weight: edgeWeight,
          };
          setEdges((prev) => [...prev, newEdge]);
        }
        setPendingFrom(null);
      }
    } else if (mode === "selectStart") {
      setStartNode(nodeId);
    }
  };

  /**
   * Run Dijkstra from the startNode we set in "selectStart" mode.
   * Infinity means: that node is not reachable from the start node in an undirected path.
   */
  const handleRunDijkstra = () => {
    if (startNode == null) {
      alert("Please select a start node first ('Select Start' mode, then click a node).");
      return;
    }
    const adjList = buildAdjacencyList(nodes, edges);
    const result = dijkstra(adjList, startNode);
    setDistances(result);
  };

  // UI
  return (
    <div style={{ display: "flex", width: "100vw", height: "100vh" }}>
      {/* --- Sidebar --- */}
      <div style={{ width: 250, borderRight: "1px solid #ccc", padding: 10, overflowY: "auto" }}>
        <h2>Dijkstra Demo</h2>
        <p>
          <strong>Mode:</strong> {mode}
        </p>
        <div style={{ marginBottom: 8 }}>
          <button onClick={() => { setMode("addNode"); setPendingFrom(null); }}>Add Node</button>{" "}
          <button onClick={() => { setMode("addEdge"); setPendingFrom(null); }}>Add Edge</button>{" "}
          <button onClick={() => { setMode("selectStart"); setPendingFrom(null); }}>Select Start</button>
        </div>

        <div style={{ marginBottom: 12 }}>
          <label>Edge Weight: </label>
          <input
            type="number"
            value={edgeWeight}
            onChange={(e: ChangeEvent<HTMLInputElement>) => setEdgeWeight(Number(e.target.value))}
            style={{ width: 60 }}
          />
        </div>

        <div style={{ marginBottom: 12 }}>
          <button onClick={handleRunDijkstra}>Run Dijkstra</button>
        </div>

        <div style={{ marginBottom: 12 }}>
          <p>Start Node: {startNode ?? "None"}</p>
          <p>Pending "From" Node: {pendingFrom ?? "None"}</p>
        </div>

        <div>
          <h4>Distances</h4>
          {Object.keys(distances).length === 0 && <p>(No results yet)</p>}
          <ul>
            {Object.entries(distances).map(([nodeId, dist]) => (
              <li key={nodeId}>
                Node {nodeId}: {dist === Infinity ? "∞ (unreachable)" : dist}
              </li>
            ))}
          </ul>
        </div>
      </div>

      {/* --- Main SVG for nodes/edges --- */}
      <div style={{ flexGrow: 1 }}>
        <svg
          width="100%"
          height="100%"
          style={{ background: "#fefefe" }}
          onClick={handleSvgClick}
        >
          {/* Edges */}
          {edges.map((edge) => {
            const fromN = nodes.find((n) => n.id === edge.from);
            const toN = nodes.find((n) => n.id === edge.to);
            if (!fromN || !toN) return null;
            const x1 = fromN.x;
            const y1 = fromN.y;
            const x2 = toN.x;
            const y2 = toN.y;
            const midX = (x1 + x2) / 2;
            const midY = (y1 + y2) / 2;

            return (
              <g key={edge.id} style={{ pointerEvents: "none" }}>
                <line
                  x1={x1}
                  y1={y1}
                  x2={x2}
                  y2={y2}
                  stroke="black"
                  strokeWidth={2}
                />
                <text
                  x={midX}
                  y={midY}
                  fill="red"
                  fontSize={12}
                  textAnchor="middle"
                  dominantBaseline="middle"
                >
                  {edge.weight}
                </text>
              </g>
            );
          })}

          {/* Nodes (circles) on top */}
          {nodes.map((node) => {
            const isStart = node.id === startNode;
            const isPendingFrom = node.id === pendingFrom;
            return (
              <circle
                key={node.id}
                cx={node.x}
                cy={node.y}
                r={15}
                fill={isStart ? "orange" : isPendingFrom ? "lightgreen" : "lightblue"}
                stroke="#333"
                strokeWidth={1}
                onClick={(e) => handleNodeClick(node.id, e)}
              />
            );
          })}

          {/* Node labels */}
          {nodes.map((node) => (
            <text
              key={`label-${node.id}`}
              x={node.x}
              y={node.y + 4}
              fontSize={12}
              fill="black"
              textAnchor="middle"
              pointerEvents="none"
            >
              {node.id}
            </text>
          ))}
        </svg>
      </div>
    </div>
  );
}
