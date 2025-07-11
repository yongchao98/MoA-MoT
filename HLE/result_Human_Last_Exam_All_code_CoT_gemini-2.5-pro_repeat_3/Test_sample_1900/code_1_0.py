def solve_complexity():
  """
  This function determines and prints the computational complexity for the two questions.
  """
  # For Question A, the complexity is constant time.
  # The number 1 is part of the final equation O(1).
  complexity_A_num = 1
  complexity_A = f"O({complexity_A_num})"

  # For Question B, the complexity is polynomial time.
  # A known algorithm runs in O(n^3). The number 3 is part of the final equation.
  complexity_B_num = 3
  complexity_B = f"O(n^{complexity_B_num})"

  # The final answer is the two complexities separated by a semicolon.
  print(f"{complexity_A}; {complexity_B}")

solve_complexity()