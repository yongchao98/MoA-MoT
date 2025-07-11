import math

def u_r(n):
  """
  Calculates the minimal order of the Picard-Fuchs differential equation
  for the given Hamiltonian system.
  The formula is based on established results for this class of potentials.
  u_r(n) = floor((n-1)/2).
  """
  return math.floor((n - 1) / 2)

def solve_and_print():
  """
  Solves for u_r(n) for n from 3 to 12 and prints the results.
  """
  results = []
  print("Calculating the values for {u_r(3), u_r(4), ..., u_r(12)}:")
  for n in range(3, 13):
    order = u_r(n)
    results.append(order)
    print(f"For n = {n}: u_r({n}) = floor(({n}-1)/2) = {order}")
  
  print("\nThe requested set is:")
  print(f"{{{', '.join(map(str, results))}}}")

solve_and_print()