import sympy
from sympy import S, log, Rational

def format_frac(f):
    """Helper function to format sympy fractions for printing."""
    f = S(f)
    if not isinstance(f, Rational):
        f = Rational(f)
    if f.q == 1:
        return str(f.p)
    else:
        return f'{f.p}/{f.q}'

# Step 1: Define the problem variables
p_true = [S(1)/2, S(1)/4, S(1)/8, S(1)/8]
q_incorrect = [S(1)/4, S(1)/2, S(1)/8, S(1)/8]
# Payout odds are 'x-for-1', meaning a 1 unit bet returns x units.
odds = [S(4), S(3), S(7), S(7)]

# Step 2: Calculate the optimal strategy and growth rate (W*)
# Check for favorable bets using true probabilities: p_i * o_i > 1
# p1*o1 = 1/2 * 4 = 2 > 1 (Favorable)
# p2*o2 = 1/4 * 3 = 3/4 < 1 (Unfavorable)
# p3*o3 = 1/8 * 7 = 7/8 < 1 (Unfavorable)
# p4*o4 = 1/8 * 7 = 7/8 < 1 (Unfavorable)
# The optimal strategy is to bet only on Bike 1.
# The optimal fraction b1* solves d/db[p1*ln(1+(o1-1)b) + (1-p1)ln(1-b)] = 0
# For p1=1/2, o1=4: d/db[1/2*ln(1+3b) + 1/2*ln(1-b)] = 0 => b1_star = 1/3.
b_star_fractions = [S(1)/3, S(0), S(0), S(0)]
B_star_total = sum(b_star_fractions)
# Calculate wealth factors S* for the optimal strategy
S_star_factors = [(1 - B_star_total + b_star_fractions[i] * odds[i]) for i in range(4)]
# Calculate optimal growth rate W*
W_star = sum(p_true[i] * log(S_star_factors[i]) for i in range(4))


# Step 3: Determine the user's betting strategy based on incorrect probabilities (q)
# Check for favorable bets using incorrect probabilities: q_i * o_i > 1
# q1*o1 = 1/4 * 4 = 1 (Neutral)
# q2*o2 = 1/2 * 3 = 3/2 > 1 (Favorable)
# q3*o3 = 1/8 * 7 = 7/8 < 1 (Unfavorable)
# q4*o4 = 1/8 * 7 = 7/8 < 1 (Unfavorable)
# The user will bet only on Bike 2.
# The optimal fraction b2 for the user solves d/db[q2*ln(1+(o2-1)b) + (1-q2)ln(1-b)] = 0
# For q2=1/2, o2=3: d/db[1/2*ln(1+2b) + 1/2*ln(1-b)] = 0 => b2 = 1/4.
b_user_fractions = [S(0), S(1)/4, S(0), S(0)]
B_user_total = sum(b_user_fractions)

# Step 4: Calculate the actual achieved growth rate (W)
# Use user's bets (b_user_fractions) but with true probabilities (p_true)
# Calculate wealth factors S for the user's strategy
S_user_factors = [(1 - B_user_total + b_user_fractions[i] * odds[i]) for i in range(4)]
# Calculate the actual achieved growth rate W
W_achieved = sum(p_true[i] * log(S_user_factors[i]) for i in range(4))

# Step 5: Calculate the decrease in growth rate (Delta W)
delta_W = W_star - W_achieved
# Simplify the expression for delta_W to be in terms of ln(2) and ln(3)
delta_W_simplified = sympy.expand_log(delta_W, force=True).simplify()


# Step 6: Print the results in the requested format
print("The doubling rate W you will achieve is:")
p_str = [format_frac(p) for p in p_true]
s_user_str = [format_frac(s) for s in S_user_factors]
w_eq_str = f"W = ({p_str[0]}) * ln({s_user_str[0]}) + ({p_str[1]}) * ln({s_user_str[1]}) + ({p_str[2]}) * ln({s_user_str[2]}) + ({p_str[3]}) * ln({s_user_str[3]})"
print(w_eq_str)
print() # Newline for separation

# Extract coefficients of ln(2) and ln(3) for the delta_W equation
c2 = delta_W_simplified.coeff(log(2))
c3 = delta_W_simplified.coeff(log(3))

c2_str = format_frac(c2)
c3_str = format_frac(abs(c3))
sign = "-" if c3 < 0 else "+"

delta_w_eq_str = f"ΔW = W* - W = ({c2_str}) * ln(2) {sign} ({c3_str}) * ln(3)"

print("By how much has your doubling rate decreased, ΔW = W* - W?")
print(delta_w_eq_str)

final_answer_W = f"W = ({p_str[0]})*ln({s_user_str[0]}) + ({p_str[1]})*ln({s_user_str[1]}) + ({p_str[2]})*ln({s_user_str[2]}) + ({p_str[3]})*ln({s_user_str[3]})"
final_answer_Delta_W = f"ΔW = ({c2_str})*ln(2) - ({c3_str})*ln(3)"
final_answer = f"{final_answer_W}, and {final_answer_Delta_W}"
# The <<<>>> format usually expects one answer, but since the question asks for two,
# I will combine them. Let's provide the final calculated forms.
print(f'<<<{final_answer}>>>')