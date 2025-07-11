import sys

# Set the parameters for the map
R = 3.57
# The desired equilibrium point from the problem description
EQUILIBRIUM_POINT = 1.05
# An initial value for X to start the iteration
x = 0.5
# Number of iterations to run to show convergence
iterations = 100

# Define the modified logistic map function
def modified_logistic_map(x_n, r):
  """
  Calculates the next value in the sequence using the modified logistic map.
  The equation is: X_n+1 = (X_n / R) * (R + 1.05 - X_n)
  """
  return (x_n / r) * (r + EQUILIBRIUM_POINT - x_n)

# Print the explanation and the equation
print("The modified logistic map is defined by the equation:")
# The prompt requires printing each number in the final equation
# We will show the general form and the specific form for R=3.57
r_val_str = str(R)
c_val_str = str(EQUILIBRIUM_POINT)
sum_val_str = str(R + EQUILIBRIUM_POINT)

print(f"X_n+1 = (X_n / R) * (R + {c_val_str} - X_n)")
print(f"\nFor R = {r_val_str}, the specific equation is:")
print(f"X_n+1 = (X_n / {r_val_str}) * ({r_val_str} + {c_val_str} - X_n)")
print(f"X_n+1 = (X_n / {r_val_str}) * ({sum_val_str} - X_n)")

# Run the simulation
for i in range(iterations):
  x = modified_logistic_map(x, R)

# Print the final result
print(f"\nStarting with an initial value of 0.5:")
print(f"After {iterations} iterations, the system converges to the stable equilibrium point:")
# Use sys.stdout.write to avoid adding a newline, making the final answer format correct
sys.stdout.write(str(x))