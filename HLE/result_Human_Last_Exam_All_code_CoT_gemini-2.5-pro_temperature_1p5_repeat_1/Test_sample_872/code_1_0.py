import math

def gcd(a, b):
    """Helper function to compute the greatest common divisor."""
    return math.gcd(a, b)

def print_fraction(num, den, label=""):
    """Helper function to print labeled fractions."""
    print(f"{label}: {num} / {den}")

print("Calculating the maximum win probability by starting in the center.\n")

# The computer has two types of responses to our center opening:
# playing in a corner (4/8 prob) or an edge (4/8 prob).
prob_o_corner_num, prob_o_corner_den = 1, 2
prob_o_edge_num, prob_o_edge_den = 1, 2

print("Step 1: Calculate the win probability if the computer plays in a corner.")
# If the computer plays a corner, our optimal response leads to a situation where:
# - With probability 2/6, the computer's next move lets us create a guaranteed fork and win.
# - With probability 2/6, the computer's next move forces us to block, after which we create a threat.
#   The computer has a 1/4 chance to block. If it fails (3/4 chance), we win.
#   If it blocks, we create another threat. The computer has a 1/2 chance to block. If it fails (1/2 chance), we win.
#   This path's win probability is 3/4 + (1/4 * 1/2) = 7/8.
# - With probability 2/6, the computer's next move again lets us create a guaranteed fork and win.
# Total win probability for this case = (2/6)*1 + (2/6)*(7/8) + (2/6)*1 = 2/3 + 7/24 = 23/24.
p_win_if_o_corner_num, p_win_if_o_corner_den = 23, 24
print_fraction(p_win_if_o_corner_num, p_win_if_o_corner_den, "P(Win | Computer plays corner)")

print("\nStep 2: Calculate the win probability if the computer plays on an edge.")
# If the computer plays an edge, our optimal response creates a threat.
# - The computer must block (1/6 chance). If it fails (5/6 chance), we win.
# - If it blocks, we create another threat. The computer must block (1/4 chance). If it fails (3/4 chance), we win.
# - If it blocks again, we create a final threat. The computer must block (1/2 chance). If it fails (1/2 chance), we win.
# - If it blocks all three times, the game is a draw (a loss for us).
# Total win probability for this case = 5/6 + (1/6 * 3/4) + (1/6 * 1/4 * 1/2) = 47/48.
p_win_if_o_edge_num, p_win_if_o_edge_den = 47, 48
print_fraction(p_win_if_o_edge_num, p_win_if_o_edge_den, "P(Win | Computer plays on edge)")

print("\nStep 3: Combine the probabilities for the final result.")
# Total P(Win) = P(O plays corner) * P(Win | O plays corner) + P(O plays edge) * P(Win | O plays edge)
#              = (1/2) * (23/24) + (1/2) * (47/48)
#              = 23/48 + 47/96
contrib_corner_num = prob_o_corner_num * p_win_if_o_corner_num
contrib_corner_den = prob_o_corner_den * p_win_if_o_corner_den
contrib_edge_num = prob_o_edge_num * p_win_if_o_edge_num
contrib_edge_den = prob_o_edge_den * p_win_if_o_edge_den

# Summing the fractions: 23/48 + 47/96 = 46/96 + 47/96
total_num = (contrib_corner_num * 2) + contrib_edge_num
total_den = contrib_edge_den

print(f"\nThe final equation is ({contrib_corner_num}/{contrib_corner_den}) + ({contrib_edge_num}/{contrib_edge_den})")
print(f"= ({contrib_corner_num * 2}/{total_den}) + ({contrib_edge_num}/{total_den})")
print(f"= {total_num}/{total_den}")

# Simplify the fraction
common_divisor = gcd(total_num, total_den)
final_num = total_num // common_divisor
final_den = total_den // common_divisor

print(f"\nAfter reducing the fraction, the maximum chance of winning is:")
print(f"{final_num} / {final_den}")

# Output the required format
# print(f"<<<{final_num}/{final_den}>>>")