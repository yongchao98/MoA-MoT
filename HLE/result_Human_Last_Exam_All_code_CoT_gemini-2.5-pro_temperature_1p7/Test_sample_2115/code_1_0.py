import math

# Step 1: Define the integral based on the assumption that u(x,y,-y,1) = u(x,y,-y,0).
# The integral to solve is: ∫[from 0 to 1] -3 * (2*e^(2x) + e^x) / (e^(2x) + e^x + 1) dx
# This integral is of the form -3 * ∫ f'(x)/f(x) dx, where f(x) = e^(2x) + e^x + 1.
# The solution is -3 * [ln(f(x))] from 0 to 1.
# Which evaluates to: -3 * (ln(f(1)) - ln(f(0))) = 3 * (ln(f(0)) - ln(f(1)))

print("Step-by-step calculation of the integral:")
print("The analytical solution to the integral is 3 * (ln(e^0 + e^0 + 1) - ln(e^2 + e + 1))")
print("This simplifies to 3 * (ln(3) - ln(e^2 + e + 1)), or 3 * ln(3 / (e^2 + e + 1)).\n")

# Step 2: Calculate the values of the terms in the expression.
e = math.exp(1)
val_at_0 = math.pow(e, 2*0) + math.pow(e, 0) + 1
val_at_1 = math.pow(e, 2*1) + math.pow(e, 1) + 1

print(f"Let f(x) = e^(2x) + e^x + 1.")
print(f"Value of e = {e:.5f}")
print(f"Evaluating at the lower limit, x=0: f(0) = e^0 + e^0 + 1 = {val_at_0}")
print(f"Evaluating at the upper limit, x=1: f(1) = e^2 + e + 1 = {val_at_1:.5f}\n")

# Step 3: Calculate the components of the final formula.
ln_val_at_0 = math.log(val_at_0)
ln_val_at_1 = math.log(val_at_1)
print(f"The integral is equal to 3 * (ln(f(0)) - ln(f(1))).")
print(f"ln(f(0)) = ln({val_at_0}) = {ln_val_at_0:.5f}")
print(f"ln(f(1)) = ln({val_at_1:.5f}) = {ln_val_at_1:.5f}\n")


# Step 4: Compute and print the final result.
final_result = 3 * (ln_val_at_0 - ln_val_at_1)
print("Final Equation:")
print(f"Result = 3 * ({ln_val_at_0:.5f} - {ln_val_at_1:.5f})")
print(f"Result = {final_result:.7f}")
