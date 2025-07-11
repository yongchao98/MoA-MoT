import math

# Define the constants from the derived formula
term1_numerator = 7
pi = math.pi
term1_subtrahend = 25
denominator_I = 6

# I_1 and I_2 calculation results
I1_val_part1 = -10/3
I1_val_part2_num = 7 * pi
I1_val_part2_den = 6
I1 = I1_val_part1 + I1_val_part2_num / I1_val_part2_den
print(f"The integral of T1(sqrt(2)x) from 0 to 1 is I_1 = -10/3 + 7*pi/6 = {I1}")

I2_min = -5/6
print(f"The minimum integral of T2(x) from 0 to 1 is I_2 = -5/6 = {I2_min}")

# Calculate the total minimum energy E_min
# E_min = 1 + 0.5 * (I_1 + I_2)
E_min_numerator = term1_numerator * pi - 13
E_min_denominator = 12

E_min = E_min_numerator / E_min_denominator

print(f"\nThe minimum total energy is given by the formula: (7 * pi - 13) / 12")
print(f"E_min = (7 * {pi} - 13) / 12")
print(f"E_min = ({term1_numerator * pi} - 13) / 12")
print(f"E_min = ({E_min_numerator}) / 12")
print(f"Final calculated value: {E_min}")
<<<0.7492552326703964>>>