import math

def modified_logistic_map(x, r):
  """
  This is the modified logistic map.
  The standard map is X_n+1 = R * X_n * (1 - X_n).
  This modified version is X_n+1 = R * X_n / (1 + X_n)^2.
  """
  # Avoid division by zero if x is -1, though in this context x is positive.
  if 1 + x == 0:
    return float('inf')
  return r * x / ((1 + x)**2)

def solve_and_demonstrate():
  """
  Solves for the equilibrium point and demonstrates the map's convergence.
  """
  R = 3.57
  
  # Calculate the theoretical equilibrium point
  # X* = sqrt(R) - 1
  equilibrium_point = math.sqrt(R) - 1
  
  # Print the modified equation
  print("The modified logistic map equation is: X_n+1 = R * X_n / (1 + X_n)^2")
  print(f"For R = {R}, the standard map is chaotic.")
  print("With the modified map, we get a stable equilibrium point.")
  
  # Print the calculated equilibrium point
  print("\nThe calculated equilibrium point X* is sqrt(R) - 1:")
  print(f"X* = sqrt({R}) - 1 = {equilibrium_point}")
  
  # Demonstrate convergence
  print("\nDemonstrating convergence from a starting point X_0 = 0.5:")
  x_n = 0.5
  print(f"X_0 = {x_n}")
  for i in range(1, 11):
    x_n_plus_1 = modified_logistic_map(x_n, R)
    print(f"X_{i} = {x_n_plus_1}")
    x_n = x_n_plus_1
    
  print(f"\nThe map converges to the equilibrium point of approximately {x_n:.6f}.")

solve_and_demonstrate()

# The final answer is the value of the equilibrium point.
final_answer = math.sqrt(3.57) - 1
# The final output still needs to print the numbers in the equation
# X* = sqrt(3.57) - 1
print(f"\nFinal Equation: {final_answer:.4f} = sqrt({3.57}) - {1}")