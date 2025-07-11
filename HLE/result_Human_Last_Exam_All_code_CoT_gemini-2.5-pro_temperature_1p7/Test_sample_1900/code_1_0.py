import math

def solve():
  """
  This function prints the computational complexity for the two questions.
  """

  # For Question A: What is the computational complexity of deciding whether I can draw a line through all the red balls in the game?
  # The problem states properties of the set of red balls. Let's analyze these properties.
  # 1. The set of red balls forms a connected component. In graph terms, the graph is connected.
  # 2. For any ball (vertex) A, the subgraph induced by its neighbors N_A is connected.
  # A graph with property 2 must be 2-connected (it cannot have a single articulation point).
  # A theorem by Shi Ronghua (1989) states that a 2-connected graph where the neighborhood of every vertex is connected is Hamiltonian.
  # This means a Hamiltonian cycle, and thus a Hamiltonian path (a line through all balls), is guaranteed to exist.
  # Since the answer to the decision question is always "yes", the complexity is constant time.
  # We represent this as O(1).
  complexity_A_constant = 1
  
  # For Question B: If I know that I can draw such a line, what is the computational complexity of finding a line that goes through all the red balls?
  # This is the search problem for the Hamiltonian Path.
  # While the existence of a path is guaranteed, the proof is non-constructive, meaning it doesn't directly provide an efficient algorithm to find the path.
  # The problem of finding a Hamiltonian Path is NP-hard in general.
  # Without a known specialized polynomial-time algorithm for this specific class of graphs, the complexity is determined by the best general-purpose algorithms.
  # The Held-Karp algorithm solves this problem using dynamic programming with a time complexity of O(n^2 * 2^n).
  # This represents an exponential complexity. The numbers in the final equation are 2 and 2.
  complexity_B_base = 2
  complexity_B_power = 2

  print(f"O({complexity_A_constant}); O(n^{complexity_B_power} * {complexity_B_base}^n)")

solve()
<<<O(1); O(n^2 * 2^n)>>>