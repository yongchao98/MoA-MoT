from fractions import Fraction

# Step 1: Determine n_0 and m_0 from analysis.
# Analysis of the problem data (as described in the plan) reveals the following:
# - The Q values are: Q_7(12.7W) > Q_1(11.6W) > Q_5(8.4W) > Q_3(7.6W) > Q_9(7.2W).
# - This corresponds to materials: Plots 1 & 7 are Copper, Plots 3 & 5 are Platinum, and Plot 9 is Carbon Steel.
# - Thus, the carbon steel fin is in plot n_0 = 9.
# - Comparing temperature drops (a proxy for the parameter 'm') shows that square fins have steeper drops than circular fins.
# - The plot for n_0 = 9 has the steepest temperature drop, consistent with it being made of the lowest conductivity material (Carbon Steel) and having the geometry that gives a larger 'm' (square).
# - Therefore, the fin is square, and m_0 = -1.

n_0 = 9
m_0 = -1

# Step 2: Calculate R(c) and R(s)
# The general ratio of heat transfer rates is R = Q_conv / Q_adi.
# The formula is R = (1/tanh(mL)) * (tanh(mL) + h/(mk)) / (1 + (h/(mk))tanh(mL)).
# For both circular and square geometries, the given conditions lead to h/(mk) = 1.
# This simplifies the ratio to R = 1 / tanh(mL).

def get_R_from_val(val):
    """Calculates R = 1 / tanh(ln(val))"""
    # tanh(ln(x)) = (x^2 - 1) / (x^2 + 1)
    tanh_val = Fraction(val**2 - 1, val**2 + 1)
    return 1 / tanh_val

# For the circular fin, mL = ln(13)
val_c = 13
R_c = get_R_from_val(val_c)

# For the square fin, mL = ln(2)
val_s = 2
R_s = get_R_from_val(val_s)

# Step 3: Compute the final value
# The expression to compute is n_0 * (R(c)/R(s))^m_0
ratio_R = R_c / R_s
final_value = n_0 * (ratio_R) ** m_0

# Step 4: Print the components and the final result as an equation.
print(f"From problem analysis, we determine n_0 = {n_0} and m_0 = {m_0}.")
print(f"For the circular fin, R(c) = 1 / tanh(ln({val_c})) = {R_c.numerator}/{R_c.denominator}")
print(f"For the square fin, R(s) = 1 / tanh(ln({val_s})) = {R_s.numerator}/{R_s.denominator}")
print("\nFinal calculation:")
print(f"n_0 * (R(c)/R(s))^m_0 = {n_0} * (({R_c.numerator}/{R_c.denominator}) / ({R_s.numerator}/{R_s.denominator}))^({m_0})")
print(f"= {n_0} * ({ratio_R.numerator}/{ratio_R.denominator})^({m_0})")
print(f"= {n_0} * {Fraction(1/ratio_R)}")
print(f"= {final_value.numerator}/{final_value.denominator}")