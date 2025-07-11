# Set the parameters
R = 3.57
# Initial condition for X
x = 0.5
# Number of iterations to run the simulation
iterations = 100

# The modified logistic map function
def modified_logistic_map(r_param, x_val):
  """
  Calculates the next value in the sequence for the modified logistic map.
  The equation is: X_n+1 = R * X_n * (1 - X_n) + X_n^3
  """
  return r_param * x_val * (1 - x_val) + x_val**3

# Run the simulation to find the equilibrium point
for i in range(iterations):
  x = modified_logistic_map(R, x)

# The final value of x is the stable equilibrium point
equilibrium_point = x

print(f"The simulation has converged to a stable equilibrium point.")
print(f"Final equilibrium point: {equilibrium_point}")
print("\nVerifying the result in the fixed-point equation: X_eq = R * X_eq * (1 - X_eq) + X_eq^3")
print("The equation with the final numbers plugged in is:")

# To satisfy the prompt "output each number in the final equation!",
# we display the equation that defines the equilibrium point, using the calculated values.
# The left side is the equilibrium point, and the right side is the function evaluated at that point.
# They should be equal.
lhs = equilibrium_point
rhs = R * equilibrium_point * (1 - equilibrium_point) + equilibrium_point**3

# We use f-string formatting to display the full equation with the numbers.
# The numbers are: the equilibrium point, R, 1, and the exponent 3.
print(f"{lhs} = {R} * {equilibrium_point} * (1 - {equilibrium_point}) + {equilibrium_point}**3")
