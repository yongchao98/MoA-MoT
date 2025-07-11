from fractions import Fraction

# Step 1: Define n_0 and m_0 based on the analysis of the problem.
# n_0 is the plot number for the carbon steel fin. Our analysis identifies this as plot 3.
n_0 = 3
# m_0 is -1 if the fin geometry is square, and 1 if circular. Plot 3 is a square fin.
m_0 = -1

# Step 2: Calculate R(c) and R(s) using the provided conditions.
# The formula for the ratio R simplifies to 1/tanh(mL) under the given conditions.
# We use the identity tanh(ln(z)) = (z^2 - 1) / (z^2 + 1).

# For the circular fin, R(c) = 1 / tanh(ln(13))
tanh_ln13_num = 13**2 - 1
tanh_ln13_den = 13**2 + 1
tanh_ln13_frac = Fraction(tanh_ln13_num, tanh_ln13_den)
Rc_frac = 1 / tanh_ln13_frac

# For the square fin, R(s) = 1 / tanh(ln(2))
tanh_ln2_num = 2**2 - 1
tanh_ln2_den = 2**2 + 1
tanh_ln2_frac = Fraction(tanh_ln2_num, tanh_ln2_den)
Rs_frac = 1 / tanh_ln2_frac

# Step 3: Compute the final expression.
ratio_frac = Rc_frac / Rs_frac
final_value_frac = n_0 * (ratio_frac ** m_0)

# Output the numbers and the final calculation as requested.
print("Final Calculation Steps:")
print(f"n_0 = {n_0}")
print(f"m_0 = {m_0}")
print(f"R(c) = 1 / tanh(ln(13)) = 1 / ({tanh_ln13_num}/{tanh_ln13_den}) = {Rc_frac.numerator}/{Rc_frac.denominator}")
print(f"R(s) = 1 / tanh(ln(2)) = 1 / ({tanh_ln2_num}/{tanh_ln2_den}) = {Rs_frac.numerator}/{Rs_frac.denominator}")
print("\nFinal Expression: n_0 * (R(c)/R(s))^m_0")
print(f"= {n_0} * (({Rc_frac.numerator}/{Rc_frac.denominator}) / ({Rs_frac.numerator}/{Rs_frac.denominator}))^({m_0})")
print(f"= {n_0} * ({ratio_frac.numerator}/{ratio_frac.denominator})^({m_0})")
print(f"= {n_0} * ({ratio_frac.denominator}/{ratio_frac.numerator})")
print(f"= {n_0 * ratio_frac.denominator}/{ratio_frac.numerator}")
print(f"The final value is: {final_value_frac.numerator}/{final_value_frac.denominator}")
