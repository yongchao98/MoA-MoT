from fractions import Fraction

# Step 1: Define the sum of the poles from the B(z) equation.
# The poles are the roots of 4z^4 - z^3 + z^2 + 1 = 0.
# By Vieta's formulas, the sum of the roots is -(-1)/4.
a4 = 4
a3 = -1
sum_poles_B = -Fraction(a3, a4)
num_poles_B = 4

# Step 2: Define the poles from the E(z) equation.
# The poles are located at z=1 and z=2.
poles_E = [Fraction(1), Fraction(2)]
sum_poles_E = sum(poles_E)
num_poles_E = len(poles_E)

# Step 3: Calculate the total sum and total number of poles.
total_sum_poles = sum_poles_B + sum_poles_E
total_num_poles = num_poles_B + num_poles_E

# Step 4: Calculate the average value.
average_z = total_sum_poles / total_num_poles

# Print the calculation steps
print(f"The sum of the 4 poles from B(z) is: {sum_poles_B.numerator}/{sum_poles_B.denominator}")
print(f"The poles from E(z) are at z=1 and z=2. Their sum is: {sum_poles_E.numerator}/{sum_poles_E.denominator}")
print(f"The total sum of all {total_num_poles} poles is: {sum_poles_B.numerator}/{sum_poles_B.denominator} + {sum_poles_E.numerator}/{sum_poles_E.denominator} = {total_sum_poles.numerator}/{total_sum_poles.denominator}")
print(f"The average value of the coordinates is the total sum divided by the number of poles.")
print(f"Average z = ({total_sum_poles.numerator}/{total_sum_poles.denominator}) / {total_num_poles} = {average_z.numerator}/{average_z.denominator}")

# Print the final result in float format for the final answer block.
print(f"The numerical value is: {float(average_z)}")
