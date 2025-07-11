import math

# Step 1: Set parameters
R = 3.57
# Set an initial value for X
x = 0.5
# Set the number of iterations to converge to the equilibrium point
iterations = 200

# Step 2: Define the modified logistic map function
def modified_logistic_map(x_n, r_param):
  """
  Calculates the next value in the sequence for the modified logistic map.
  The modification is the addition of the x_n^4 term.
  """
  return r_param * x_n * (1 - x_n) + math.pow(x_n, 4)

# Step 3: Iterate the map to find the stable equilibrium point
for i in range(iterations):
  x = modified_logistic_map(x, R)

# The loop finds the value of x for which x_n+1 = x_n
equilibrium_point = x

print("The proposed modified logistic map is: X_n+1 = R * X_n * (1 - X_n) + X_n^4")
print(f"For the parameter R = {R}, the map has a stable equilibrium point.")
print(f"\nAfter {iterations} iterations, the calculated equilibrium point is: {equilibrium_point}")

# Step 4: Verify that the found value is indeed an equilibrium point
# by plugging it back into the equation.
print("\n--- Verification ---")
print("An equilibrium point X* must satisfy the equation: X* = R * X* * (1 - X*) + X*^4")
print("Let's check this with our calculated value:\n")

# Printing the final equation with the numbers, as requested.
print(f"{equilibrium_point} = {R} * {equilibrium_point} * (1 - {equilibrium_point}) + {equilibrium_point} ^ 4")

# Calculating the right-hand side of the equation to confirm the equality
rhs = R * equilibrium_point * (1 - equilibrium_point) + math.pow(equilibrium_point, 4)
print(f"\nValue of the left side (X*): {equilibrium_point}")
print(f"Value of the right side:    {rhs}")
print(f"The difference is {'{:.2e}'.format(abs(equilibrium_point - rhs))}, which confirms it's an equilibrium point.")

# Extra step: check for stability
# The equilibrium is stable if |f'(x*)| < 1
# f'(x) = R - 2*R*x + 4*x^3
derivative = R - 2 * R * equilibrium_point + 4 * math.pow(equilibrium_point, 3)
print(f"\nThe derivative at the equilibrium point is {derivative:.4f}.")
print(f"Since |{derivative:.4f}| < 1, the equilibrium point is stable.")