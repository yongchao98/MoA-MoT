import math

# Step 1: Define the constants lambda1 and lambda2 from the problem description.
# lambda1 = (1 + sqrt(17))/2
# lambda2 = (1 - sqrt(17))/2
lambda1 = (1 + math.sqrt(17)) / 2
lambda2 = (1 - math.sqrt(17)) / 2

# Step 2: Set the initial condition x2(0).
# The problem is framed as a "controllability problem". In this context, if the initial
# state is not provided, it's standard to assume the system starts from the origin, x(0) = 0.
# Therefore, we set x2(0) = 0.
x2_0 = 0

# The original expression is:
# ((2/3)*lambda1*exp(lambda1/2) - (1/3)*exp(lambda1/2))*x2(0) - (2/3)*lambda2*exp(lambda2/2) - (10/3)*exp(lambda1/2)
# With x2(0) = 0, the first term is zero, and the expression simplifies.

# Step 3: Calculate the numerical value of each component in the simplified expression.
# The simplified expression is: - (2/3) * lambda2 * exp(lambda2/2) - (10/3) * exp(lambda1/2)
c1 = -2/3
val_lambda2 = lambda2
val_exp_lambda2_half = math.exp(lambda2 / 2)
c2 = -10/3
val_exp_lambda1_half = math.exp(lambda1 / 2)

# Step 4: Print the final equation with the numbers substituted, as requested.
print("Based on the standard assumption that x(0)=0, the expression simplifies.")
print("The final equation is of the form: (C1 * V1 * V2) + (C2 * V3)")
print(f"C1 = {c1}")
print(f"V1 (lambda2) = {val_lambda2}")
print(f"V2 (exp(lambda2/2)) = {val_exp_lambda2_half}")
print(f"C2 = {c2}")
print(f"V3 (exp(lambda1/2)) = {val_exp_lambda1_half}")
print("\nSubstituting these values into the equation:")
print(f"Value = ({c1}) * ({val_lambda2}) * ({val_exp_lambda2_half}) + ({c2}) * ({val_exp_lambda1_half})")

# Step 5: Calculate and print the final result.
final_value = (c1 * val_lambda2 * val_exp_lambda2_half) + (c2 * val_exp_lambda1_half)
print(f"\nFinal calculated value: {final_value}")

print(f"\n<<<{final_value}>>>")