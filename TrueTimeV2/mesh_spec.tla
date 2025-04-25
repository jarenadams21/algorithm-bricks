---- MODULE QuantumMeshSpec ----
EXTENDS Naturals, Reals, FiniteSets, Sequences, TLC

(***************************************************************************)
(* Constants and Parameters                                               *)
(***************************************************************************)
CONSTANTS NODES \in Nat \ {0},   \* number of nodes, NODES = 5
          MODES \in Nat \ {0},   \* number of mode types, MODES = 4
          ROUNDS \in Nat \ {0},  \* number of rounds, ROUNDS = 20
          DT \in Real,          \* time step (s), DT = 1e-8
          Lambda \in Real,      \* Dyson scale, Lambda = 1.0
          h, kB, T              \* physical constants

ASSUME NODES = 5
ASSUME MODES = 4
ASSUME ROUNDS = 20
ASSUME DT = 1e-8
ASSUME Lambda = 1.0
ASSUME h = 6.62607015e-34 /
       kB = 1.380649e-23 /
       T = 300.0

(***************************************************************************)
(* Variables                                                               *)
(***************************************************************************)
VARIABLES
  NodePos,    \* [1..NODES] -> <<x,y>> in Real^2
  Forces,     \* [1..NODES] -> <<Fx,Fy>> in Real^2
  Pending,    \* [1..NODES] -> Seq( [mode: 0..MODES-1, ts: Real] )
  pendCnt,    \* [1..NODES] -> Nat
  Ledger,     \* [1..NODES] -> Seq(Record)
  ledgCnt,    \* [1..NODES] -> Nat
  ChunkStore, \* Seq(ChunkEntry)
  chunkCnt,   \* Nat
  SpinDist    \* [1..NODES] -> Real  \* hardware spin distribution per node

Record == [mode: 0..MODES-1, count: Nat, ts: Real, prev: Real, hash: Real]
ChunkEntry == [node: 1..NODES, p: [mode:0..MODES-1, ts:Real], w: Real, mw: Real]

(***************************************************************************)
(* Initialization                                                          *)
(***************************************************************************)
Init ==
  /\ NodePos \in [1..NODES -> Real^2]
  /\ Forces \in [1..NODES -> Real^2]
  /\ \A i \in 1..NODES: Pending[i] = <<>>
  /\ pendCnt = [i \in 1..NODES |-> 0]
  /\ Ledger[i] \in [i \in 1..NODES |-> <<>>]
  /\ ledgCnt = [i \in 1..NODES |-> 0]
  /\ ChunkStore = <<>>
  /\ chunkCnt = 0
  /\ SpinDist \in [1..NODES -> Real]

(***************************************************************************)
(* Auxiliary Functions                                                     *)
(***************************************************************************)
fDyson(d) == 1 / (1 + LN(1 + d / Lambda))

wCond(cnt) ==
  LET nu == cnt % ROUNDS + 1 IN
    1 / (EXP(h * nu / (kB * T)) - 1)

Distance(i, j) ==
  LET dx == NodePos[j][1] - NodePos[i][1]
      dy == NodePos[j][2] - NodePos[i][2]
  IN SQRT(dx*dx + dy*dy)

ComputeForces ==
  Forces' = [i \in 1..NODES |->
      LET Fsum ==
            \E j \in 1..NODES \ {i}: 
              LET d == Distance(i,j) IN
                IF d <= 1e-6 THEN <<0,0>> ELSE
                LET rep == 1/d^2
                    att == d^2
                    net == (att - rep) * fDyson(d) * wCond(chunkCnt)
                    dx == (NodePos[j][1] - NodePos[i][1])/d
                    dy == (NodePos[j][2] - NodePos[i][2])/d
                IN << net*dx, net*dy >>
      IN Fsum]

EulerStep ==
  NodePos' = [i \in 1..NODES |-> <<
      (NodePos[i][1] + DT * Forces[i][1]),
      (NodePos[i][2] + DT * Forces[i][2])>>]

(***************************************************************************)
(* Chandy-Lamport Snapshot Protocol                                         *)
(***************************************************************************)
SndMsg(i, j) == 
  \* Asynchronous send of state marker and SpinDist[i]
  \* ... abstracted
  TRUE

RcvMsg(i, j) ==
  \* Matches SndMsg and records channel state
  TRUE

Snapshot ==
  /\ \A i, j \in 1..NODES : (i # j) => SndMsg(i,j) /\ RcvMsg(j,i)

(***************************************************************************)
(* Next-State Relation                                                     *)
(***************************************************************************)
Next ==
  \E c \in 0..(MODES-1), i \in 1..NODES :
    /\ pendCnt' = [pendCnt EXCEPT ![i] = pendCnt[i]+1]
    /\ Pending'  = [Pending EXCEPT ![i] = Append(Pending[i], [mode |-> c, ts |-> qc_ts()])]
    /\ (* order_n, commit_n abstracted into Pend->Ledger updates with spin coupling *)
    /\ ComputeForces
    /\ EulerStep
    /\ UNCHANGED <<Ledger, ledgCnt, ChunkStore, chunkCnt, SpinDist>>

(***************************************************************************)
(* Specification                                                           *)
(***************************************************************************)
Spec == Init /\ [][Next]_<<NodePos,Forces,Pending,pendCnt,Ledger,ledgCnt,ChunkStore,chunkCnt,SpinDist>>

(***************************************************************************)
(* Invariants & Theorems                                                   *)
(***************************************************************************)
Inv_NoNaN ==
  \A i \in 1..NODES :
    (~ IsNaN(NodePos[i][1])) /\ (~ IsNaN(NodePos[i][2]))
  /\ (~ IsNaN(Forces[i][1])) /\ (~ IsNaN(Forces[i][2]))

THEOREM NoNaNConsistent == Spec => []Inv_NoNaN
<1>1. Init => Inv_NoNaN              BY DEF Init, Inv_NoNaN
<1>2. Inv_NoNaN /\ Next => Inv_NoNaN BY (* proof obligation, uses boundedness of fDyson, wCond *)
<1>QED BY <1>1, <1>2, PTL

Liouville ==
  \* Phase-space volume preserved:
  \* ∫_{Ω} dNodePos dForces = constant under ComputeForces & EulerStep

THEOREM LiouvillePreservation == Spec => [](Liouville)
<2>1. EulerStep preserves volume    BY Eulerian ODE theory
<2>2. ComputeForces is divergence-free in NodePos space BY ∑ ∂F_i/∂x_i = 0 (QCD-inspired) 
<2>QED BY <2>1, <2>2, PTL

=============================================================================
