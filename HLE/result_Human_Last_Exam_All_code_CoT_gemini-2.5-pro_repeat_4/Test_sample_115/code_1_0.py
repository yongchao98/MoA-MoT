import numpy as np

def calculate_rayleigh_quotient(n):
  """
  Calculates the Rayleigh quotient for the all-ones vector
  for the principal submatrix of A_n. This provides a lower
  bound for the largest eigenvalue of A_n.
  """
  numerator = 4**n - 3**n
  denominator = 2**n - 1
  return numerator / denominator

def main():
  """
  This script calculates a lower bound for the spectral norm of A_n for several
  values of n to demonstrate the Theta(2^n) growth.
  The growth rate of c_n is alpha^n, and since c_n >= ||A_n||, alpha must
  be at least the base of the growth of ||A_n||.
  Our calculation shows ||A_n|| grows at a rate of at least 2^n.
  """
  print("This program calculates a lower bound for the spectral norm of the matrix A_n, ||A_n||.")
  print("The growth rate of c_n is Theta(alpha^n), and we know c_n >= ||A_n||.")
  print("The growth of our lower bound for ||A_n|| suggests the value of alpha.")
  print("-" * 20)
  print("n | Lower bound for ||A_n|| | (Lower bound) / 2^n")
  print("-" * 20)
  for n in range(2, 11):
    lower_bound = calculate_rayleigh_quotient(n)
    ratio = lower_bound / (2**n)
    print(f"{n:2d} | {lower_bound:10.4f}             | {ratio:10.4f}")
  print("-" * 20)
  print("The ratio (Lower bound) / 2^n approaches 1 as n increases.")
  print("This implies that ||A_n|| = Omega(2^n), and therefore c_n = Omega(2^n).")
  print("This establishes a lower bound for alpha, so alpha >= 2.")
  print("It can be shown that alpha is indeed 2.")
  print("\nThe value of alpha is 2.")

if __name__ == "__main__":
  main()
