from fractions import Fraction

# Step 1: Define the determined parameters n_0 and m_0.
# Based on analysis of the problem description and image data:
# n_0 is the plot number for the carbon steel fin, which is 9.
# m_0 specifies the geometry of the carbon steel fin (1 for circular, -1 for square), which is circular.
n_0 = 9
m_0 = 1

# Step 2: Calculate R(c) for the circular fin.
# Given conditions lead to R(c) = 1 / tanh(ln(13)).
# tanh(ln(x)) = (x^2 - 1) / (x^2 + 1)
tanh_ln13 = Fraction(13**2 - 1, 13**2 + 1)
R_c = 1 / tanh_ln13

# Step 3: Calculate R(s) for the square fin.
# Given conditions lead to R(s) = 1 / tanh(ln(2)).
tanh_ln2 = Fraction(2**2 - 1, 2**2 + 1)
R_s = 1 / tanh_ln2

# Step 4: Calculate the final expression: n_0 * (R(c) / R(s))^m_0
ratio_R = R_c / R_s
final_value = n_0 * (ratio_R ** m_0)

# Step 5: Print the final equation with all numbers, and the result in x/y format.
print(f"The calculation is performed as follows:")
print(f"{n_0} * (({R_c.numerator}/{R_c.denominator}) / ({R_s.numerator}/{R_s.denominator}))^{m_0} = {final_value.numerator}/{final_value.denominator}")