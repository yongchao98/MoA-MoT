from fractions import Fraction

# Step 1: Identify n_0 and m_0 from image analysis.
# The heat transfer rate 'q' depends on thermal conductivity 'k'. For the materials given, k_Copper > k_Platinum > k_CarbonSteel.
# The heat transfer values from the plots are: q7=12.74, q1=11.63, q5=8.43, q3=7.63, q9=7.19 (in W).
# The lowest value, q=7.19 W from plot 9, must correspond to the Carbon Steel fin. Therefore, n_0 = 9.
n_0 = 9
# To determine the geometry, we match the tip temperature of plot 9 (just under 50 °C) with the 3D plots.
# Plot 4 shows a square fin with a tip temperature of 48.6596 °C, which is a clear match.
# The geometry is square, so m_0 is -1.
m_0 = -1

# Step 2: Calculate R(c), the ratio for the circular fin.
# The general formula is R = (1 + (h/mk) * coth(mL)) / (1 + (h/mk) * tanh(mL)).
# For the circular fin, we are given conditions that simplify to mL = ln(13) and h/mk = 1.
# We calculate tanh(ln(13)) = (13^2 - 1) / (13^2 + 1) = 168 / 170.
# So, tanh(ln(13)) = Fraction(168, 170) which simplifies to 84/85.
tanh_ln13 = Fraction(84, 85)
# coth(ln(13)) is the reciprocal.
coth_ln13 = 1 / tanh_ln13
# R(c) = (1 + coth(ln(13))) / (1 + tanh(ln(13)))
Rc = (1 + coth_ln13) / (1 + tanh_ln13)

# Step 3: Calculate R(s), the ratio for the square fin.
# For the square fin, the given conditions simplify to mL = ln(2) and h/mk = 1.
# We calculate tanh(ln(2)) = (2^2 - 1) / (2^2 + 1) = 3 / 5.
tanh_ln2 = Fraction(3, 5)
# coth(ln(2)) is the reciprocal.
coth_ln2 = 1 / tanh_ln2
# R(s) = (1 + coth(ln(2))) / (1 + tanh(ln(2)))
Rs = (1 + coth_ln2) / (1 + tanh_ln2)

# Step 4: Compute the final value n_0 * (R(c)/R(s))^m_0.
ratio_Rc_Rs = Rc / Rs
final_term = ratio_Rc_Rs ** m_0
result = n_0 * final_term

# Step 5: Print the final answer with the full equation showing each number.
print(f"The expression to be evaluated is: n_0 * (R(c)/R(s))^m_0")
print("Substituting the determined and calculated values, the equation becomes:")
print(f"{n_0} * (({Rc.numerator}/{Rc.denominator}) / ({Rs.numerator}/{Rs.denominator}))^({m_0})")
print(f"= {n_0} * ({ratio_Rc_Rs.numerator}/{ratio_Rc_Rs.denominator})^({m_0})")
print(f"= {n_0} * ({final_term.numerator}/{final_term.denominator})")
print(f"= {result.numerator}/{result.denominator}")