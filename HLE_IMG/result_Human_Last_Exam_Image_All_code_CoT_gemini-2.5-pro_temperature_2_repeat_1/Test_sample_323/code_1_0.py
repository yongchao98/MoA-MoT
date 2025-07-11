from fractions import Fraction

# Step 1: Define n_0 and m_0 based on our analysis.
# n_0 is the plot number for the carbon steel fin. We identified this as plot 9.
n_0 = 9
# m_0 is 1 for a circular fin and -1 for a square fin. Plot 9 is a circular fin.
m_0 = 1

# Step 2: Calculate R(c) for the circular fin.
# The conditions given for the circular fin are L/(k/h) = 4L/d = ln(13).
# From this, we can show that for this specific problem, h/(m*k) = 1.
# The parameter mL = ln(13).
# The ratio R simplifies to R = 1 / tanh(mL) when h/(m*k)=1.
# We use the identity tanh(ln(x)) = (x^2 - 1) / (x^2 + 1).
x_c = 13
tanh_ml_c = Fraction(x_c**2 - 1, x_c**2 + 1)  # (169-1)/(169+1) = 168/170 = 84/85
Rc = 1 / tanh_ml_c  # 1 / (84/85) = 85/84

# Step 3: Calculate R(s) for the square fin.
# The conditions given for the square fin are L/(k/h) = 4L/w = ln(2).
# Similarly, for this case h/(m*k) = 1 and mL = ln(2).
x_s = 2
tanh_ml_s = Fraction(x_s**2 - 1, x_s**2 + 1)  # (4-1)/(4+1) = 3/5
Rs = 1 / tanh_ml_s # 1 / (3/5) = 5/3

# Step 4: Compute the final expression.
# The expression is n_0 * (R(c)/R(s))^m_0
final_result = n_0 * (Rc / Rs)**m_0

# Print the final equation with all values substituted
print("Final Equation:")
print(f"{n_0} * (({Rc.numerator}/{Rc.denominator}) / ({Rs.numerator}/{Rs.denominator}))^{m_0}")

# Print the final answer in the format x/y
print("\nFinal Answer:")
print(f"{final_result.numerator}/{final_result.denominator}")
