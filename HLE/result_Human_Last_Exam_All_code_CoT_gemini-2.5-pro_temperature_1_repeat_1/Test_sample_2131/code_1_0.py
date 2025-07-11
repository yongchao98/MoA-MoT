# The analytical solution for the deflection y at x=0 is y(0) = (3/4) * (9/4)^(3/5).
# This script calculates the numerical value of this expression.

# Define the numbers in the final equation
numerator_1, denominator_1 = 3, 4
numerator_2, denominator_2 = 9, 4
numerator_3, denominator_3 = 3, 5

# Perform the calculation step-by-step
fraction_1 = numerator_1 / denominator_1
base = numerator_2 / denominator_2
exponent = numerator_3 / denominator_3
power_result = base ** exponent
final_result = fraction_1 * power_result

# Output the equation and the steps of the calculation
print(f"The final equation for the deflection y(0) is: y(0) = ({numerator_1}/{denominator_1}) * ({numerator_2}/{denominator_2})^({numerator_3}/{denominator_3})")
print("\n--- Calculation Steps ---")
print(f"1. Calculate the base of the power: {numerator_2}/{denominator_2} = {base}")
print(f"2. Calculate the exponent: {numerator_3}/{denominator_3} = {exponent}")
print(f"3. Evaluate the power: ({base})^({exponent}) = {power_result}")
print(f"4. Perform the final multiplication: {fraction_1} * {power_result} = {final_result}")
print("\n-------------------------")
print(f"\nThe membrane's deflection at x = 0, y(0), is: {final_result}")
