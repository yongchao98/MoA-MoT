import math

def modified_logistic_map(r, x):
  """
  The modified logistic map function.
  X_n+1 = R * X_n * (1 - X_n) + 0.3 * R * X_n^3
  """
  term1 = r * x * (1 - x)
  term2 = 0.3 * r * (x**3)
  return term1 + term2

# Set the parameter R
R = 3.57
# Set the initial value for X
x0 = 1.0

print(f"Investigating the modified logistic map for R = {R}")
print(f"The equation is: X_n+1 = R * X_n * (1 - X_n) + 0.3 * R * X_n^3")
print(f"Starting with X_0 = {x0}\n")

x = x0
# Iterate to demonstrate convergence to the new equilibrium point
for i in range(20):
    x_next = modified_logistic_map(R, x)
    print(f"X_{i+1:2d} = {x_next:.6f}")
    # Stop if converged
    if abs(x_next - x) < 1e-9:
        x = x_next
        break
    x = x_next

print(f"\nSystem converges to a stable equilibrium point at approximately: {x:.6f}\n")

# Per the instructions, output each number in the final equation.
# The general form is: X_n+1 = r_val * X_n * (one_val - X_n) + c_val * r_val * X_n^k_val
r_val = 3.57
one_val = 1
c_val = 0.3
k_val = 3

print("Final equation with numerical values:")
print(f"X_n+1 = {r_val} * X_n * ({one_val} - X_n) + {c_val} * {r_val} * X_n^{k_val}")

# To calculate the precise value for the final answer
# we solve the fixed point equation: x = f(x)
# x = R*x*(1-x) + 0.3*R*x^3
# For x != 0: 1 = R*(1-x) + 0.3*R*x^2
# 0.3*R*x^2 - R*x + (R-1) = 0
a = 0.3 * R
b = -R
c = R - 1
# Quadratic formula: x = (-b +/- sqrt(b^2 - 4ac)) / 2a
discriminant = b**2 - 4 * a * c
sol1 = (-b + math.sqrt(discriminant)) / (2 * a)
sol2 = (-b - math.sqrt(discriminant)) / (2 * a)
# We choose the solution close to 1.05
final_answer = sol2
print(f"\n# The precise equilibrium point is calculated as {final_answer:.6f}")
<<<1.051012>>>