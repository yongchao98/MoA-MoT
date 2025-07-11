import numpy as np

def h_function(x):
  """
  A sample smooth function h: R^2 -> R that is in L^1.
  Here h(x1, x2) = exp(-(x1^2 + x2^2)).
  """
  return np.exp(-(x[0]**2 + x[1]**2))

def get_kronecker_delta(i, j):
  """
  Computes the Kronecker delta. Returns 1 if i==j, 0 otherwise.
  Indices are 1-based.
  """
  return 1 if i == j else 0

def calculate_question_mark_1(h_func, x, i, j):
  """
  Calculates the term ?_1 based on the derived formula.
  ?_1 = (1/2) * h(x) * delta_ij
  """
  h_val = h_func(x)
  delta_ij = get_kronecker_delta(i, j)
  term_value = 0.5 * h_val * delta_ij
  
  # The task asks to output each number in the final equation.
  print(f"For i={i}, j={j}:")
  print(f"?_1 = (1/2) * h(x) * delta_ij")
  print(f"?_1 = {0.5} * {h_val:.4f} * {delta_ij}")
  print(f"?_1 = {term_value:.4f}")
  print("-" * 20)

# Define a point x
x_point = (0.5, 1.0)

print(f"Calculating ?_1 for h(x) = exp(-(x1^2 + x2^2)) at x = {x_point}\n")

# Loop through all combinations of i and j (1 or 2)
for i_idx in [1, 2]:
  for j_idx in [1, 2]:
    calculate_question_mark_1(h_function, x_point, i_idx, j_idx)
