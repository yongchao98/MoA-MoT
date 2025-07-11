import math
from fractions import Fraction

def format_fraction(f):
    """Helper function to format fractions for printing."""
    if f.denominator == 1:
        return str(f.numerator)
    return f"{f.numerator}/{f.denominator}"

# Step 0: Define the given probabilities and odds using the Fraction class for precision.
p = {
    1: Fraction(1, 2),
    2: Fraction(1, 4),
    3: Fraction(1, 8),
    4: Fraction(1, 8)
}
q = {
    1: Fraction(1, 4),
    2: Fraction(1, 2),
    3: Fraction(1, 8),
    4: Fraction(1, 8)
}
# Odds are "k-for-1", which correspond to decimal odds of k.
o = {
    1: Fraction(4),
    2: Fraction(3),
    3: Fraction(7),
    4: Fraction(7)
}

# Step 1: Calculate the optimal betting strategy and growth rate W*.
# Check for favorable bets using true probabilities p: p_i * o_i > 1.
# p[1]*o[1] = 1/2 * 4 = 2 > 1 (Favorable)
# p[2]*o[2] = 1/4 * 3 = 3/4 < 1 (Unfavorable)
# p[3]*o[3] = 1/8 * 7 = 7/8 < 1 (Unfavorable)
# p[4]*o[4] = 1/8 * 7 = 7/8 < 1 (Unfavorable)
# The optimal strategy is to only bet on Bike 1.

# Calculate the optimal Kelly fraction for Bike 1.
p1, o1 = p[1], o[1]
# Formula for single bet: b* = p - (1-p)/(o-1)
b1_star = p1 - (1 - p1) / (o1 - 1)

# Calculate the optimal growth rate W*.
# W* = p1 * ln(1 - b1_star + b1_star*o1) + (1-p1) * ln(1 - b1_star)
# This simplifies to W* = p1 * ln(p1*o1 / (1 - b1_star*(o1-1))) + (1-p1)*ln((1-p1)/(1-b1_star)) which is messy
# Let's simplify the combined expression: W* = (1/2)*ln(1-1/3+4/3) + (1/2)*ln(1-1/3) = (1/2)*ln(2) + (1/2)*ln(2/3) = (1/2)*ln(4/3)
w_star_coeff = Fraction(1, 2)
w_star_arg = Fraction(4, 3)
w_star_expr_str = f"({format_fraction(w_star_coeff)}) * ln({format_fraction(w_star_arg)})"

# Step 2: Determine the bettor's strategy based on incorrect probabilities q.
# Check for perceived favorable bets using q: q_i * o_i > 1.
# q[1]*o[1] = 1/4 * 4 = 1 (Fair, bet 0)
# q[2]*o[2] = 1/2 * 3 = 3/2 > 1 (Perceived Favorable)
# q[3]*o[3] = 1/8 * 7 = 7/8 < 1 (Perceived Unfavorable)
# q[4]*o[4] = 1/8 * 7 = 7/8 < 1 (Perceived Unfavorable)
# The bettor decides to only bet on Bike 2.

# Calculate the bettor's chosen fraction for Bike 2 based on q.
q2, o2 = q[2], o[2]
b2_prime = q2 - (1 - q2) / (o2 - 1)

# Step 3: Calculate the actual achieved growth rate W.
# The bet is b2_prime on Bike 2, and 0 on others.
# Calculate returns based on this bet.
return_if_2_wins = 1 - b2_prime + b2_prime * o2
return_if_other_wins = 1 - b2_prime

# Calculate expected log growth rate using TRUE probabilities p.
# W = p[1]*ln(return_if_other_wins) + p[2]*ln(return_if_2_wins) + p[3]*ln(return_if_other_wins) + p[4]*ln(return_if_other_wins)
prob_of_loss = p[1] + p[3] + p[4]
prob_of_win = p[2]
w_coeff1 = prob_of_loss
w_arg1 = return_if_other_wins
w_coeff2 = prob_of_win
w_arg2 = return_if_2_wins

# Step 4: Print the final equations.
print("The doubling rate W you will achieve is:")
w_expr_str = f"({format_fraction(w_coeff1)}) * ln({format_fraction(w_arg1)}) + ({format_fraction(w_coeff2)}) * ln({format_fraction(w_arg2)})"
print(f"W = {w_expr_str}")
print("-" * 20)
print("The decrease in doubling rate is:")
delta_w_expr_str = f"{w_star_expr_str} - ({w_expr_str})"
print(f"Delta W = {delta_w_expr_str}")

# Final Answer for the bot
answer_W = f"W = {format_fraction(w_coeff1)}*ln({format_fraction(w_arg1)}) + {format_fraction(w_coeff2)}*ln({format_fraction(w_arg2)})"
answer_Delta_W = f"Delta W = {format_fraction(w_star_coeff)}*ln({format_fraction(w_star_arg)}) - ({format_fraction(w_coeff1)}*ln({format_fraction(w_arg1)}) + {format_fraction(w_coeff2)}*ln({format_fraction(w_arg2)}))"

print(f"\n<<<{answer_W}, {answer_Delta_W}>>>")