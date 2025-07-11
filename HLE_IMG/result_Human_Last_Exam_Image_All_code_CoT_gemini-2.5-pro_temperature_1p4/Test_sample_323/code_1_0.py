from fractions import Fraction

# 1. Set the determined parameters
n_0 = 9
m_0 = -1

# 2. Define a function to calculate R using fractions for precision
def calculate_R(a, b):
    """
    Calculates R = (1/tanh(mL)) * (tanh(mL) + h/(mk)) / (1 + (h/(mk)) * tanh(mL))
    where mL = ln(a) and h/(mk) = b.
    """
    # tanh(ln(a)) = (a^2 - 1) / (a^2 + 1)
    tanh_mL = Fraction(a**2 - 1, a**2 + 1)
    
    numerator = tanh_mL + b
    denominator = 1 + b * tanh_mL
    
    R = (1 / tanh_mL) * (numerator / denominator)
    return R

# 3. Calculate R(c) for the circular fin
# Given conditions lead to a=13 and b=1/2
a_c = 13
b_val = Fraction(1, 2)
R_c = calculate_R(a_c, b_val)

# 4. Calculate R(s) for the square fin
# Given conditions lead to a=2 and b=1/2
a_s = 2
R_s = calculate_R(a_s, b_val)

# 5. Calculate the final expression: n_0 * (R_c / R_s)^m_0
# This simplifies to n_0 * (R_s / R_c) since m_0 = -1
final_value = Fraction(n_0) * (R_s / R_c)

# 6. Print the results in the required format
print(f"Based on the problem analysis, we have:")
print(f"n_0 = {n_0}")
print(f"m_0 = {m_0}")
print(f"\nUnder the specified conditions, the calculated ratios are:")
print(f"R(c) = {R_c.numerator}/{R_c.denominator}")
print(f"R(s) = {R_s.numerator}/{R_s.denominator}")
print(f"\nThe value of the expression n_0 * (R(c)/R(s))^m_0 is calculated as:")
print(f"{n_0} * ({R_c.numerator}/{R_c.denominator} / {R_s.numerator}/{R_s.denominator}) ^ {m_0}")
print(f"= {n_0} * ({R_s.numerator}/{R_s.denominator} / {R_c.numerator}/{R_c.denominator})")
print(f"= {n_0} * ({R_s.numerator * R_c.denominator} / {R_s.denominator * R_c.numerator})")
print(f"= {final_value.numerator}/{final_value.denominator}")