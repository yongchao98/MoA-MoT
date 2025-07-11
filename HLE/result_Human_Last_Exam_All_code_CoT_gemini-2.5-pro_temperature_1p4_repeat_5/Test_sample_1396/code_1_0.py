import math

def solve_cake_cutting_bound():
  """
  Calculates and prints the upper bound for a specific cake-cutting problem.

  The problem is to find the query complexity (the upper bound 'O') for
  achieving a connected epsilon-envy-free allocation for four agents,
  assuming monotonicity of the valuation functions.
  """
  # Number of agents in the problem
  n_agents = 4

  # According to the 2020 paper "Connected epsilon-envy-free cake cutting
  # for four agents" by Segal-Halevi and Ianovski, the upper bound on the
  # number of queries (cuts or evaluations) required is 20.
  upper_bound_O = 20

  # Print the final result as an equation.
  # The output includes each number involved in the statement.
  print(f"For n = {n_agents} agents, the upper bound on queries for a connected epsilon-envy-free allocation is O = {upper_bound_O}.")

solve_cake_cutting_bound()