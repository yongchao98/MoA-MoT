import math
from fractions import Fraction

# Step 1: Determine n_0 and m_0 from the problem description and image analysis.
# Based on thermal properties, the carbon steel fin corresponds to the plot with the
# lowest heat transfer rate and lowest tip temperature for a given geometry.
# Comparing square fins (plots 1, 3, 5), plot 3 has the lowest q and T_L.
n_0 = 3
# The geometry for plot 3 is square.
m_0 = -1
print(f"Step 1: Identification of parameters")
print(f"The carbon steel fin is in plot n_0 = {n_0}")
print(f"The geometry is square, so m_0 = {m_0}\n")

# Step 2: Calculate R(c) for the circular fin
# For the circular fin, mL = ln(13) and h/(mk) = 1.
# The ratio R = Q_conv / Q_adi simplifies to 1/tanh(mL).
x_c = 13
tanh_ln_xc_num = x_c - 1/x_c
tanh_ln_xc_den = x_c + 1/x_c
# Using Fraction for exact rational arithmetic
tanh_ln_13 = Fraction(int(tanh_ln_xc_num * x_c), int(tanh_ln_xc_den * x_c))
Rc = 1 / tanh_ln_13

print(f"Step 2: Calculation of R(c) for the circular fin")
print(f"R(c) = 1 / tanh(ln({x_c}))")
print(f"tanh(ln({x_c})) = ({x_c} - 1/{x_c}) / ({x_c} + 1/{x_c}) = {tanh_ln_13.numerator}/{tanh_ln_13.denominator}")
print(f"R(c) = 1 / ({tanh_ln_13.numerator}/{tanh_ln_13.denominator}) = {Rc.numerator}/{Rc.denominator}\n")

# Step 3: Calculate R(s) for the square fin
# For the square fin, mL = ln(2) and h/(mk) = 1.
# The ratio R simplifies to 1/tanh(mL).
x_s = 2
tanh_ln_xs_num = x_s - 1/x_s
tanh_ln_xs_den = x_s + 1/x_s
# Using Fraction for exact rational arithmetic
tanh_ln_2 = Fraction(int(tanh_ln_xs_num * x_s), int(tanh_ln_xs_den * x_s))
Rs = 1 / tanh_ln_2

print(f"Step 3: Calculation of R(s) for the square fin")
print(f"R(s) = 1 / tanh(ln({x_s}))")
print(f"tanh(ln({x_s})) = ({x_s} - 1/{x_s}) / ({x_s} + 1/{x_s}) = {tanh_ln_2.numerator}/{tanh_ln_2.denominator}")
print(f"R(s) = 1 / ({tanh_ln_2.numerator}/{tanh_ln_2.denominator}) = {Rs.numerator}/{Rs.denominator}\n")

# Step 4: Calculate the final expression
# The expression is n_0 * (R(c)/R(s))^m_0
ratio_R = Rc / Rs
final_result = n_0 * (ratio_R)**m_0

print(f"Step 4: Final Calculation")
print(f"The expression to evaluate is: n_0 * (R(c)/R(s))^m_0")
print(f"Substituting values: {n_0} * (({Rc.numerator}/{Rc.denominator}) / ({Rs.numerator}/{Rs.denominator}))^({m_0})")
print(f"= {n_0} * ({ratio_R.numerator}/{ratio_R.denominator})^({m_0})")
print(f"= {n_0} * {Fraction(ratio_R.denominator, ratio_R.numerator)}")
print(f"= {final_result.numerator}/{final_result.denominator}")

print(f"\nFinal Answer:")
print(f"{final_result.numerator}/{final_result.denominator}")