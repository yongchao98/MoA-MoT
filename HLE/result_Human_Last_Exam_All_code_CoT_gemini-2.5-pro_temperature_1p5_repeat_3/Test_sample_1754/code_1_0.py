import math

# Step 1: Define the given constants from the problem.
# lambda1 and lambda2 are defined as per the problem statement.
# x2(0) is determined to be 0 from our analysis of a simple valid solution.
sqrt17 = math.sqrt(17)
lambda1 = (1 + sqrt17) / 2
lambda2 = (1 - sqrt17) / 2
x2_0 = 0.0
t = 0.5  # time t = 1/2

# Step 2: Calculate the value of each term in the full expression.
# Full Expression: TermA * x2(0) + TermB + TermC

# TermA is the coefficient of x2(0)
termA = ((2/3) * lambda1 * math.exp(lambda1 * t) - (1/3) * math.exp(lambda1 * t))

# TermB is the term with lambda2
termB = - (2/3) * lambda2 * math.exp(lambda2 * t)

# TermC is the last term
termC = - (10/3) * math.exp(lambda1 * t)

# Step 3: Calculate the final value using x2(0) = 0.
final_value = termA * x2_0 + termB + termC

# As requested, output each number in the final equation.
# Since x2(0) is 0, the final equation simplifies to the sum of TermB and TermC.
print("--- Calculation Breakdown ---")
print(f"Based on a simplified equilibrium solution, we found x2(0) = {x2_0}")
print(f"The equation to be calculated is:")
print(f"((2/3 * λ1 * e^(λ1/2) - 1/3 * e^(λ1/2)) * x2(0)) - (2/3 * λ2 * e^(λ2/2)) - (10/3 * e^(λ1/2))")
print("\n--- Values of the components ---")
print(f"Value of the term multiplying x2(0): {termA}")
print(f"Value of the second term (-2/3 * λ2 * e^(λ2/2)): {termB}")
print(f"Value of the third term (-10/3 * e^(λ1/2)): {termC}")

# Step 4: Print the final result.
print("\n--- Final Result ---")
print(f"Final Value = ({termA}) * {x2_0} + ({termB}) + ({termC})")
print(f"Final Value = {final_value}")
