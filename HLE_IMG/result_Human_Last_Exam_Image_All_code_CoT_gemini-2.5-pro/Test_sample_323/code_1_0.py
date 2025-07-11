from fractions import Fraction

def calculate_tanh_ln_x(x):
    """Calculates tanh(ln(x)) and returns it as a Fraction."""
    # tanh(y) = (e^y - e^-y) / (e^y + e^-y)
    # If y = ln(x), then e^y = x and e^-y = 1/x
    # tanh(ln(x)) = (x - 1/x) / (x + 1/x)
    numerator = x - (1/x)
    denominator = x + (1/x)
    return Fraction(numerator).limit_denominator() / Fraction(denominator).limit_denominator()

# Part 1: Determined from problem analysis
n_0 = 9
m_0 = 1

# Part 2: Calculate R(c) and R(s)
# For the circular fin case
val_c = 13
tanh_ln_c = calculate_tanh_ln_x(val_c)
R_c = 1 / tanh_ln_c

# For the square fin case
val_s = 2
tanh_ln_s = calculate_tanh_ln_x(val_s)
R_s = 1 / tanh_ln_s

# Part 3: Final Calculation
final_value = n_0 * (R_c / R_s)**m_0

# Print the final equation with all numerical values
# The format is n_0 * (R(c) / R(s))^m_0 = result
print("The final equation is:")
print(f"{n_0} * (({R_c.numerator}/{R_c.denominator}) / ({R_s.numerator}/{R_s.denominator}))^{m_0} = {final_value.numerator}/{final_value.denominator}")

# Print the final answer in the requested format
print("\nThe final answer is:")
print(f"<<<{final_value.numerator}/{final_value.denominator}>>>")
