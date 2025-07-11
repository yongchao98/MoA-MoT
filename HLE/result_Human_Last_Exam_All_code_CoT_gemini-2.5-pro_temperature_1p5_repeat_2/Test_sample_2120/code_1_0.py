from fractions import Fraction

# Poles from B(z)
# The sum of roots for 4z^4 - z^3 + z^2 + 1 = 0 is -(-1)/4 = 1/4. Number of roots is 4.
# The sum of roots for z^4 + z^2 - z + 4 = 0 is -0/1 = 0. Number of roots is 4.
sum_poles_B = Fraction(1, 4)
num_poles_B = 8

# Poles from E(z)
# From the analysis of the determinant, the sum of poles is -c5/c6.
# c6 = 550, c5 = -4900. Sum = 4900/550 = 98/11. Number of poles is 6.
sum_poles_E = Fraction(98, 11)
num_poles_E = 6

# Total sum and number of poles
total_sum_of_poles = sum_poles_B + sum_poles_E
total_number_of_poles = num_poles_B + num_poles_E

# Calculate the average
average_value = total_sum_of_poles / total_number_of_poles

# Print the final equation
print("The average value is calculated by summing all pole coordinates and dividing by the number of poles.")
print("Average = (Sum of B's poles + Sum of E's poles) / (Number of B's poles + Number of E's poles)")
print(f"Average = (({sum_poles_B.numerator}/{sum_poles_B.denominator}) + ({sum_poles_E.numerator}/{sum_poles_E.denominator})) / ({num_poles_B} + {num_poles_E})")
print(f"Average = (({total_sum_of_poles.numerator}/{total_sum_of_poles.denominator})) / {total_number_of_poles}")
print(f"Average = {average_value.numerator}/{average_value.denominator}")

# Final Answer
# print(f'<<<{average_value.numerator}/{average_value.denominator}>>>')