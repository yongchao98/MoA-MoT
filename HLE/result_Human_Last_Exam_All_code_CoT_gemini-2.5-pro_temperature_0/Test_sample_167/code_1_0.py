def calculate_alon_tarsi_number_k_n_n(n):
  """
  Calculates the Alon-Tarsi number for a complete bipartite graph K_n,n.

  The Alon-Tarsi number for K_n,n is n + 1.
  """
  result = n + 1
  # The problem asks to output each number in the final equation.
  print(f"The Alon-Tarsi number of K_{n},{n} is calculated as: {n} + 1 = {result}")

# For the graph K_1000,1000, the value of n is 1000.
n_val = 1000
calculate_alon_tarsi_number_k_n_n(n_val)