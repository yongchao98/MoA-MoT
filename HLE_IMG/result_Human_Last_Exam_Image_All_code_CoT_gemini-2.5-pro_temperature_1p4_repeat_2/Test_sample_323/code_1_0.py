from fractions import Fraction

# Step 1: Assign values for n_0 and m_0 based on analysis of the problem and plots.
# n_0 is the plot number for the Carbon Steel fin. Based on thermal conductivity,
# Carbon Steel has the lowest heat transfer, which corresponds to Plot 9.
n_0 = 9
# m_0 depends on the geometry of fin n_0. Based on analysis of the Q_conv/Q_adi ratio,
# fin 9 has the characteristics of a square fin.
m_0 = -1

# Step 2: Calculate R(c) for the circular fin under the specified conditions.
# For the circular case, it was derived that R(c) = 1 / tanh(ln(13)).
val_c = 13
# Use Fraction for exact arithmetic. tanh(x) = (e^x-e^-x)/(e^x+e^-x)
# tanh(ln(13)) = (13 - 1/13)/(13 + 1/13) = (169-1)/(169+1) = 168/170 = 84/85
tanh_c = Fraction(val_c**2 - 1, val_c**2 + 1)
R_c = 1 / tanh_c

# Step 3: Calculate R(s) for the square fin under the specified conditions.
# For the square case, it was derived that R(s) = 1 / tanh(ln(2)).
val_s = 2
# tanh(ln(2)) = (2 - 1/2)/(2 + 1/2) = (4-1)/(4+1) = 3/5
tanh_s = Fraction(val_s**2 - 1, val_s**2 + 1)
R_s = 1 / tanh_s

# Step 4: Compute the final result using the derived values.
# The expression is n_0 * (R(c)/R(s))^m_0.
ratio_Rc_Rs = R_c / R_s
final_result = n_0 * (ratio_Rc_Rs ** m_0)

# Output the components of the equation and the final answer.
print("The final equation is: n_0 * (R(c) / R(s)) ^ m_0")
print(f"n_0 = {n_0}")
print(f"m_0 = {m_0}")
print(f"R(c) = 1 / tanh(ln({val_c})) = {R_c.numerator}/{R_c.denominator}")
print(f"R(s) = 1 / tanh(ln({val_s})) = {R_s.numerator}/{R_s.denominator}")
print(f"The ratio R(c)/R(s) = {ratio_Rc_Rs.numerator}/{ratio_Rc_Rs.denominator}")
print(f"Final calculation: {n_0} * ({ratio_Rc_Rs.numerator}/{ratio_Rc_Rs.denominator}) ^ ({m_0})")
print(f"Result: {final_result.numerator}/{final_result.denominator}")