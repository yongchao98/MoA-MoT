import math

# Step 1: Define the derived function f(z) = sqrt(pi) / Gamma(z)
def f(z):
    """
    Calculates the value of the derived function f(z) = sqrt(pi) / Gamma(z).
    The Gamma function has poles at non-positive integers. Its reciprocal is 0 there.
    """
    if z <= 0 and z == int(z):
        return 0.0
    return math.sqrt(math.pi) / math.gamma(z)

# Step 2: Choose a test value for z.
# A non-integer value is chosen to avoid the zeros of 1/Gamma(z).
z_test_value = 4.5

# Step 3: Verify the functional equation: f(z) = 2^(1-z) * f(z/2) * f((z+1)/2)
print("--- Verification of the Functional Equation ---")
print(f"Testing for z = {z_test_value}")

# Calculate the left-hand side (LHS)
lhs = f(z_test_value)

# Calculate the right-hand side (RHS)
# Note: In Python, ** is the power operator.
factor = 2**(1 - z_test_value)
f_z_div_2 = f(z_test_value / 2)
f_z_plus_1_div_2 = f((z_test_value + 1) / 2)
rhs = factor * f_z_div_2 * f_z_plus_1_div_2

print(f"LHS = f(z) = {lhs}")
print(f"RHS = 2^(1-z) * f(z/2) * f((z+1)/2) = {rhs}")

# We use math.isclose to handle potential floating-point inaccuracies
if math.isclose(lhs, rhs):
    print("Result: The derived function satisfies the functional equation for the test value.")
else:
    print("Result: The derived function does not satisfy the functional equation for the test value.")

# Step 4: Verify the initial condition f(1) = sqrt(pi)
print("\n--- Verification of the Initial Condition ---")
f_1 = f(1)
sqrt_pi = math.sqrt(math.pi)
print(f"Calculated f(1) = {f_1}")
print(f"Expected value sqrt(pi) = {sqrt_pi}")
if math.isclose(f_1, sqrt_pi):
    print("Result: The derived function satisfies the initial condition f(1) = sqrt(pi).")
else:
    print("Result: The derived function does not satisfy the initial condition f(1) = sqrt(pi).")

# Step 5: Print the final explicit form of the function, outputting its numbers
print("\n--- Explicit Form of f(z) ---")
# The equation is f(z) = sqrt(pi) / Gamma(z).
# The numbers in the equation are pi and the power 0.5 for the square root.
numerator_base = math.pi
numerator_power = 0.5
print(f"The explicit form of the function is:")
# Print the equation with its numerical components
print(f"f(z) = ({numerator_base}**{numerator_power}) / Gamma(z)")
print("Or, more commonly written as: f(z) = sqrt(pi) / Γ(z)")