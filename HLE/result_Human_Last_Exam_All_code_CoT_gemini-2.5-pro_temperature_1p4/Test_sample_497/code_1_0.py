import math

print("This script solves for the capacitance 'x' that makes the equivalent capacitance of a ladder circuit independent of the number of cells.")
print("-" * 80)

print("\nStep 1: Set up the equation for the characteristic capacitance.")
print("Let 'c' be the capacitance of the capacitors in each cell.")
print("Let 'x' be the characteristic capacitance of the ladder, which is the required value for the terminating capacitor.")
print("The input capacitance of a single cell loaded with a capacitance 'x' must be equal to 'x'.")

print("\nThe input capacitance `C_in` of a single cell loaded with `C_L` is given by the formula:")
print("C_in = (c^2 + c*C_L) / (3c + 2*C_L)")

print("\nSetting C_in = C_L = x, we get the equation:")
print("x = (c^2 + c*x) / (3c + 2*x)")

print("\nRearranging this into a standard quadratic form (a*x^2 + b*x + d = 0) gives:")
print("x * (3c + 2x) = c^2 + c*x")
print("3*c*x + 2*x^2 = c^2 + c*x")
print("2*x^2 + 2*c*x - c^2 = 0")
print("-" * 80)

print("\nStep 2: Solve the quadratic equation.")
print("The final equation to solve for x is: 2*x^2 + 2*c*x - c^2 = 0")
print("We can identify the coefficients for the quadratic formula (in terms of the variable x):")
# As requested, here are the numbers/expressions for each term in the final equation.
a = 2
b_expr = "2*c"
d_expr = "-c^2"
print(f"Coefficient of x^2: a = {a}")
print(f"Coefficient of x:   b = {b_expr}")
print(f"Constant term:      d = {d_expr}")

print("\nUsing the quadratic formula x = [-b +/- sqrt(b^2 - 4ad)] / (2a):")
print("x = [-(2c) +/- sqrt((2c)^2 - 4*(2)*(-c^2))] / (2*2)")
print("x = [-2c +/- sqrt(4c^2 + 8c^2)] / 4")
print("x = [-2c +/- sqrt(12c^2)] / 4")
print("x = [-2c +/- 2*sqrt(3)*c] / 4")
print("x = c * (-1 +/- sqrt(3)) / 2")

print("\nSince capacitance must be a positive value, we take the positive root.")
final_expression = "c * (sqrt(3) - 1) / 2"
print(f"The required value for x is: {final_expression}")

# Calculate the numerical value of the coefficient
coefficient = (math.sqrt(3) - 1) / 2
print(f"\nNumerically, this is x = {coefficient:.6f} * c")
print("-" * 80)