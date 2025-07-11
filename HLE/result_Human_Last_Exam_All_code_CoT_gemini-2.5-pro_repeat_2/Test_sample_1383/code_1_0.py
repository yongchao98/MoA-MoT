import math
from fractions import Fraction

# This script calculates the achieved growth rate W and the decrease in growth rate ΔW
# based on the provided probabilities and odds.
# It uses the Fraction class to maintain precision and formats the output
# to show the derivation with fractions and natural logarithms, as requested.

# Step 1: Define the problem parameters
# True probabilities
p = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8)]
# Incorrectly believed probabilities (used for betting fractions)
q = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 8), Fraction(1, 8)]
# Decimal odds (d-for-1 means decimal odds are d)
d = [Fraction(4), Fraction(3), Fraction(7), Fraction(7)]

# Helper function to format fractions for printing
def format_fraction(frac):
    """Formats a Fraction object into a string, e.g., 1/2 or 5."""
    if frac.denominator == 1:
        return str(frac.numerator)
    return f"{frac.numerator}/{frac.denominator}"

# --- Part 1: Calculate the achieved growth rate W ---

print("### Calculation of the achieved doubling rate W ###\n")

print("You bet fractions b_i = q_i based on your incorrect beliefs.")
print("The achieved growth rate W is the expected log return, calculated with the true probabilities p_i.")
print("Formula: W = p_1*log(b_1*d_1) + p_2*log(b_2*d_2) + p_3*log(b_3*d_3) + p_4*log(b_4*d_4)")
print("Substituting b_i = q_i:")
print("W = p_1*log(q_1*d_1) + p_2*log(q_2*d_2) + p_3*log(q_3*d_3) + p_4*log(q_4*d_4)\n")

# Print the equation with all numbers substituted
w_eq_initial = (f"W = {format_fraction(p[0])}*log({format_fraction(q[0])}*{format_fraction(d[0])}) + "
                f"{format_fraction(p[1])}*log({format_fraction(q[1])}*{format_fraction(d[1])}) + "
                f"{format_fraction(p[2])}*log({format_fraction(q[2])}*{format_fraction(d[2])}) + "
                f"{format_fraction(p[3])}*log({format_fraction(q[3])}*{format_fraction(d[3])})")
print("Substituting the numbers:")
print(w_eq_initial)

# Print the equation with simplified products inside log
w_eq_simplified_log = (f"W = {format_fraction(p[0])}*log({format_fraction(q[0]*d[0])}) + "
                       f"{format_fraction(p[1])}*log({format_fraction(q[1]*d[1])}) + "
                       f"{format_fraction(p[2])}*log({format_fraction(q[2]*d[2])}) + "
                       f"{format_fraction(p[3])}*log({format_fraction(q[3]*d[3])})")
print("\nSimplifying the products inside the logarithms:")
print(w_eq_simplified_log)

# Print the final simplified expression
print("\nSince log(1) = 0, the first term vanishes.")
print("The last two terms are identical and can be combined: (1/8 + 1/8) * log(7/8) = 1/4 * log(7/8).")
final_w_expr = (f"W = {format_fraction(p[1])}*log({format_fraction(q[1]*d[1])}) + "
                f"{format_fraction(p[2]+p[3])}*log({format_fraction(q[2]*d[2])})")
print(f"\nThe final expression for W is:")
print(final_w_expr)


# --- Part 2: Calculate the decrease in growth rate ΔW ---

print("\n" + "="*50)
print("\n### Calculation of the decrease in growth rate ΔW ###\n")
print("The decrease in growth rate is ΔW = W* - W, where W* is the optimal rate achieved by betting b_i = p_i.")
print("This difference is equivalent to the Kullback-Leibler divergence D_KL(p || q).\n")
print("Formula: ΔW = p_1*log(p_1/q_1) + p_2*log(p_2/q_2) + p_3*log(p_3/q_3) + p_4*log(p_4/q_4)\n")

# Print the equation with all numbers substituted
dw_eq_initial = (f"ΔW = {format_fraction(p[0])}*log({format_fraction(p[0])}/{format_fraction(q[0])}) + "
                 f"{format_fraction(p[1])}*log({format_fraction(p[1])}/{format_fraction(q[1])}) + "
                 f"{format_fraction(p[2])}*log({format_fraction(p[2])}/{format_fraction(q[2])}) + "
                 f"{format_fraction(p[3])}*log({format_fraction(p[3])}/{format_fraction(q[3])})")
print("Substituting the numbers:")
print(dw_eq_initial)

# Print the equation with simplified ratios inside log
dw_eq_simplified_log = (f"ΔW = {format_fraction(p[0])}*log({format_fraction(p[0]/q[0])}) + "
                        f"{format_fraction(p[1])}*log({format_fraction(p[1]/q[1])}) + "
                        f"{format_fraction(p[2])}*log({format_fraction(p[2]/q[2])}) + "
                        f"{format_fraction(p[3])}*log({format_fraction(p[3]/q[3])})")
print("\nSimplifying the ratios inside the logarithms:")
print(dw_eq_simplified_log)

# Print the final simplified expression
print("\nSince log(1) = 0, the last two terms vanish.")
print("Since log(1/2) = -log(2), the second term can be rewritten.")
dw_final_step = f"ΔW = {format_fraction(p[0])}*log(2) - {format_fraction(p[1])}*log(2)"
print(dw_final_step)
dw_final_step_2 = f"ΔW = ({format_fraction(p[0])} - {format_fraction(p[1])})*log(2)"
print(dw_final_step_2)

final_dw_val = p[0] - p[1]
final_dw_expr = f"ΔW = {format_fraction(final_dw_val)}*log(2)"
print(f"\nThe final expression for ΔW is:")
print(final_dw_expr)
