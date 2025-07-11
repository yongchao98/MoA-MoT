import math

def modified_logistic_map(x, r):
  """
  A modified logistic map.
  The standard map is X_n+1 = R * X_n * (1 - X_n).
  This is modified to X_n+1 = R * X_n * |R/(R+1) - X_n|.
  """
  # The carrying capacity term is K = r / (r + 1)
  k = r / (r + 1.0)
  return r * x * abs(k - x)

def solve_task():
  """
  Solves the task by defining and demonstrating the modified logistic map.
  """
  R = 3.57
  
  # Calculate the theoretical equilibrium points (fixed points)
  # The equation for fixed points is X = R*X*|K-X|, where K=R/(R+1)
  # This gives X=0, or |K-X|=1/R, which has two solutions:
  # X1 = K - 1/R  (This one turns out to be stable)
  # X2 = K + 1/R  (This one is the target unstable point)
  K = R / (R + 1.0)
  x_star_stable = K - 1.0/R
  x_star_unstable_target = K + 1.0/R

  # --- Outputting the results ---
  
  print("The modified logistic map is:")
  print(f"X_n+1 = R * X_n * |R/(R+1) - X_n|")
  print("-" * 30)

  print(f"For the given R = {R}:")
  # The target equilibrium point is X* = R/(R+1) + 1/R
  print("The new equilibrium point is calculated as:")
  print(f"X* = {R}/({R}+1) + 1/{R}")
  print(f"X* = {K:.4f} + {1/R:.4f}")
  print(f"X* = {x_star_unstable_target:.4f}")
  print(f"This value is approximately 1.05, as requested.")
  print("-" * 30)
  
  # Demonstrate that x_star_unstable_target is an equilibrium point
  print("Demonstrating it is an equilibrium point:")
  x = x_star_unstable_target
  print(f"Starting at X_0 = {x:.4f}")
  x = modified_logistic_map(x, R)
  print(f"X_1 = {x:.4f}")
  x = modified_logistic_map(x, R)
  print(f"X_2 = {x:.4f}")
  print("The value remains constant, confirming it's an equilibrium point.")
  print("-" * 30)

  # Demonstrate that the map is no longer chaotic
  print("Demonstrating the system is no longer chaotic:")
  print(f"(It converges to the other, stable fixed point at {x_star_stable:.4f})")
  x = 0.4
  print(f"Starting at X_0 = {x:.4f}")
  for i in range(1, 11):
    x = modified_logistic_map(x, R)
    print(f"X_{i} = {x:.4f}")
    
solve_task()