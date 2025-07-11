import sys

# This problem is a theoretical computer science question about complexity.
# The code below does not solve the problem but formats the derived answer
# as requested by the user prompt.

def solve():
  """
  Solves the computational complexity questions based on graph theory principles.
  
  Question A: What is the computational complexity of deciding whether I can draw a line through all the red balls in the game?
  
  The problem can be modeled using a graph G where red balls are vertices and neighborhood connections are edges. The given properties are:
  1. G is connected.
  2. For every vertex v, the subgraph induced by its neighbors N(v) is connected. This is the definition of a locally connected graph.
  
  A locally connected graph is also necessarily claw-free. A famous theorem by Ryjáček states that any connected, locally connected, claw-free graph (with n>=3) has a Hamiltonian cycle, and therefore a Hamiltonian path.
  
  This means that under the given properties, a path through all n balls is *always* possible. The decision problem is thus trivial: the answer is always "yes", regardless of n. This requires no computation that scales with n.
  
  Complexity for A: O(1)
  
  Question B: If I know that I can draw such a line, what is the computational complexity of finding a line that goes through all the red balls?
  
  This asks for the complexity of *finding* the Hamiltonian path. The proofs that guarantee existence are often constructive, implying a polynomial-time algorithm.
  
  A common approach for 1-connected graphs is to decompose them into 2-connected components (blocks) and find paths in each. Algorithms exist for finding Hamiltonian cycles in 2-connected locally connected graphs with a complexity of O(k^3), where k is the number of vertices in the block.
  
  In the worst case, the entire graph is one large block of size n. Therefore, the overall complexity of finding the path is bounded by this case.
  
  Complexity for B: O(n^3)
  
  """
  
  # The final answer is a string representing the complexities.
  # The problem does not ask for an algorithm implementation, just the Big-O notation.
  answer_A = "O(1)"
  answer_B = "O(n^3)"
  
  print(f"{answer_A}; {answer_B}")

solve()