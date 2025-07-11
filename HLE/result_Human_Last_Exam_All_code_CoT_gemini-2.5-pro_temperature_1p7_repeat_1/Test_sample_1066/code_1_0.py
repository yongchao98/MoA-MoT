import math

def calculate_lower_bound(n, alpha):
  """
  Calculates the lower bound on the expected detection score E[S].

  Args:
    n: The number of tokens in the text.
    alpha: The average entropy per token.
  
  Returns:
    The lower bound for E[S].
  """
  pi = math.pi
  # The bound involves the constant zeta(2) = pi^2 / 6
  zeta_2 = (pi**2) / 6
  
  lower_bound = n * (alpha - math.log(zeta_2))
  
  # Print the final equation as a string representation
  print(f"The equation for the lower bound on E[S] is:")
  print(f"E[S] >= n * (alpha - ln(pi^2 / 6))")
  print("-" * 20)
  
  # Print the values used in the calculation
  print(f"Given values:")
  print(f"n (number of tokens) = {n}")
  print(f"alpha (average entropy) = {alpha:.4f}")
  print(f"pi = {pi:.6f}")
  print("-" * 20)
  
  # Print the calculated result
  print(f"Calculation:")
  print(f"E[S] >= {n} * ({alpha:.4f} - ln({pi**2 / 6:.4f}))")
  print(f"E[S] >= {n} * ({alpha:.4f} - {math.log(zeta_2):.4f})")
  print(f"E[S] >= {n} * ({alpha - math.log(zeta_2):.4f})")
  print(f"Lower bound for E[S] = {lower_bound:.4f}")
  
  return lower_bound

if __name__ == '__main__':
  # Example usage with some arbitrary but realistic values for a language model
  # Number of tokens
  n_tokens = 1000
  # Average entropy (e.g., perplexity is e^alpha, so a perplexity of ~20 gives alpha ~ 3.0)
  avg_entropy = 3.0
  
  calculate_lower_bound(n_tokens, avg_entropy)