import math
from fractions import Fraction

# Step 1: Define the true probabilities (p), the incorrect probabilities (q),
# and the decimal odds (d).
p = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8)]
q = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 8), Fraction(1, 8)]
# Odds are (4-for-1, 3-for-1, 7-for-1, 7-for-1), which correspond to decimal odds d.
d = [4, 3, 7, 7]

# Step 2: Calculate the achieved growth rate W.
# The betting fractions 'f' are based on the incorrect probabilities 'q'.
# The expected growth W is calculated using the true probabilities 'p'.
print("1. Calculation of the achieved growth rate W:")
print("W = p_1*log(q_1*d_1) + p_2*log(q_2*d_2) + p_3*log(q_3*d_3) + p_4*log(q_4*d_4)")
print(f"W = ({p[0]})*log(({q[0]})*{d[0]}) + ({p[1]})*log(({q[1]})*{d[1]}) + ({p[2]})*log(({q[2]})*{d[2]}) + ({p[3]})*log(({q[3]})*{d[3]})")
term1_val = q[0] * d[0]
term2_val = q[1] * d[1]
term3_val = q[2] * d[2]
term4_val = q[3] * d[3]
print(f"W = ({p[0]})*log({term1_val}) + ({p[1]})*log({term2_val}) + ({p[2]})*log({term3_val}) + ({p[3]})*log({term4_val})")
print("Since log(1) = 0, the first term vanishes.")
# Combine terms 3 and 4 as they are identical
print(f"W = ({p[1]})*log({term2_val}) + ({p[2]} + {p[3]})*log({term3_val})")
final_p_term_34 = p[2]+p[3]
print("The final expression for W is:")
print(f"W = ({p[1]})*log({term2_val}) + ({final_p_term_34})*log({term3_val})")


print("\n--------------------------------\n")

# Step 3: Calculate the optimal growth rate W*.
# The optimal betting fractions 'f*' are the true probabilities 'p'.
print("2. Calculation of the optimal growth rate W*:")
print("W* = p_1*log(p_1*d_1) + p_2*log(p_2*d_2) + p_3*log(p_3*d_3) + p_4*log(p_4*d_4)")
term1_star_val = p[0] * d[0]
term2_star_val = p[1] * d[1]
term3_star_val = p[2] * d[2]
print(f"W* = ({p[0]})*log({term1_star_val}) + ({p[1]})*log({term2_star_val}) + ({p[2]})*log({term3_star_val}) + ({p[3]})*log({term3_star_val})")
print("The final expression for W* is:")
print(f"W* = ({p[0]})*log({term1_star_val}) + ({p[1]})*log({term2_star_val}) + ({final_p_term_34})*log({term3_star_val})")


print("\n--------------------------------\n")

# Step 4: Calculate the decrease in growth rate ΔW.
print("3. Calculation of the decrease in growth rate ΔW = W* - W:")
print("This is also the Kullback-Leibler divergence D_KL(p || q) = sum(p_i * log(p_i / q_i)).")
print("ΔW = p_1*log(p_1/q_1) + p_2*log(p_2/q_2) + p_3*log(p_3/q_3) + p_4*log(p_4/q_4)")
ratio1 = p[0]/q[0]
ratio2 = p[1]/q[1]
ratio3 = p[2]/q[2]
ratio4 = p[3]/q[3]
print(f"ΔW = ({p[0]})*log({ratio1}) + ({p[1]})*log({ratio2}) + ({p[2]})*log({ratio3}) + ({p[3]})*log({ratio4})")
print("Since log(1) = 0, the last two terms are zero.")
print(f"ΔW = ({p[0]})*log({ratio1}) + ({p[1]})*log({ratio2})")
print(f"Using log(1/2) = -log(2):")
print(f"ΔW = ({p[0]})*log(2) - ({p[1]})*log(2)")
final_coeff = p[0] - p[1]
print(f"ΔW = ({p[0]} - {p[1]})*log(2)")
print("The final expression for ΔW is:")
print(f"ΔW = ({final_coeff})*log(2)")
