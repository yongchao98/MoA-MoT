import math

def solve_hamiltonicity_threshold():
  """
  This function calculates and prints the formula for the d-threshold for Hamiltonicity.
  """
  # The problem asks for the d-threshold for Hamiltonicity in H_n U G(n, p)
  # where delta(H_n) >= d = n/2 - eta.
  # The key parameter is k = n/2 - d = eta.
  
  # Based on the literature for randomly perturbed graphs, the threshold depends on the regime of eta.
  # For small eta, p is on the order of eta / n^2.
  # For large eta, p is on the order of (eta * ln(eta)) / n^2.
  
  # A unified formula that covers the entire range given (1/2 <= eta <= n/64) can be expressed
  # by combining these two regimes. We assume the leading constant of the threshold is 1.
  
  # The formula is p = (eta * max(1, ln(eta))) / n^2.
  # We will now print this final equation.
  
  print("The d-threshold for Hamiltonicity is given by the probability p:")
  # The numbers in the equation are 1 (implicit multiplier) and 2 (the exponent).
  # We will print them explicitly as requested.
  print("p = (eta * max(1, ln(eta))) / (n^2)")

solve_hamiltonicity_threshold()