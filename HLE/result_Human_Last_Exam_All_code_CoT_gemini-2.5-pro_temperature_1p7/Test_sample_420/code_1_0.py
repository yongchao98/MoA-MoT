import math

# The problem is solved symbolically to find the exact form.
# The side length 's' cancels out in the final ratio.
# Here we will construct and print the string representing the final answer.

# The final ratio is (pi + 6*sqrt(3) - 6) / 12.
# We need to print each number in this expression.
numerator_term_1_coeff = 1
numerator_term_1_symbol = "pi"
numerator_term_2_coeff = 6
numerator_term_2_symbol = "sqrt(3)"
numerator_term_3_coeff = -6
denominator = 12

print("The ratio of the area of D to the area of S can be written in its exact form as an equation.")
print("Ratio = (A*pi + B*sqrt(C) - D) / E")
print(f"Where the coefficients and constants are:")
print(f"A = {numerator_term_1_coeff}")
print(f"B = {numerator_term_2_coeff}")
print(f"C = {3}")
print(f"D = {abs(numerator_term_3_coeff)}")
print(f"E = {denominator}")
print("\nFinal Equation:")
print(f"Ratio = ({numerator_term_1_coeff}*pi + {numerator_term_2_coeff}*sqrt(3) - {abs(numerator_term_3_coeff)}) / {denominator}")
