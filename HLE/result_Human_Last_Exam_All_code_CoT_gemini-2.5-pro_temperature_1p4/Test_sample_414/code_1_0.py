import math
from scipy.optimize import fsolve

# --- Step 1: Define constants ---
# Half-lives in days
T_p = 12.75  # Ba-140 half-life in days
T_d = 40.28 / 24.0  # La-140 half-life in days (converted from hours)

# Decay constants (lambda) in days^-1
lambda_p = math.log(2) / T_p
lambda_d = math.log(2) / T_d

# Measured activities (signal) in kBq/mL
S0 = 1.4  # Signal at t=0 (after separation)
S14 = 2.1 # Signal at t=14 days

# --- Step 2: Solve for the relative Cherenkov efficiency k ---
# The signal S(t) = A_p(t) + k * A_d(t)
# At t=0, A_d(0) = 0, so S(0) = A_p(0). Thus, A_p_initial = 1.4
A_p_initial = S0

# At t=14 days:
# A_p(14) = A_p_initial * exp(-lambda_p * 14)
# A_d(14) = (lambda_d / (lambda_d - lambda_p)) * A_p_initial * (exp(-lambda_p * 14) - exp(-lambda_d * 14))
# S(14) = A_p(14) + k * A_d(14)
# S14 = A_p_initial * exp(-lambda_p * 14) + k * (lambda_d / (lambda_d - lambda_p)) * A_p_initial * (exp(-lambda_p * 14) - exp(-lambda_d * 14))

def solve_for_k(k):
    """Equation for the signal at t=14 days, to be solved for k."""
    term1 = A_p_initial * math.exp(-lambda_p * 14)
    term2 = k * (lambda_d / (lambda_d - lambda_p)) * A_p_initial * (math.exp(-lambda_p * 14) - math.exp(-lambda_d * 14))
    return term1 + term2 - S14

# Solve for k
k = fsolve(solve_for_k, 1.0)[0] # Initial guess for k is 1.0

# --- Step 3: Solve for the time t1 ---
# Assumption: The signal of the unseparated sample at time t1 is equal to the signal of the separated sample at 14 days.
# S_unsep(t1) = S14 = 2.1
#
# The signal of the unseparated sample is:
# S_unsep(t1) = A_p(t1) + k * A_d(t1)
# where A_p(t1) is the parent activity at time t1, which is our A_p_initial = 1.4
# and A_d(t1) is the daughter activity grown from an initially pure parent over time t1.
# A_d(t1) = (lambda_d / (lambda_d - lambda_p)) * A_p(t1) * (1 - exp(-(lambda_d-lambda_p)*t1))
#
# So, the equation to solve for t1 is:
# A_p_initial + k * (lambda_d / (lambda_d - lambda_p)) * A_p_initial * (1 - exp(-(lambda_d - lambda_p) * t1)) = S14

def solve_for_t1(t1):
    """Equation for the signal of the unseparated mixture at time t1."""
    if t1 < 0: return float('inf') # t1 must be positive
    ingrowth_factor = (lambda_d / (lambda_d - lambda_p)) * (1 - math.exp(-(lambda_d - lambda_p) * t1))
    s_unsep = A_p_initial + k * A_p_initial * ingrowth_factor
    return s_unsep - S14

# Solve for t1
t1 = fsolve(solve_for_t1, 1.0)[0] # Initial guess for t1 is 1.0 day

# --- Step 4: Output the result ---
# Let's print all the numbers in the final equation used to find t1
# The equation is S14 = A_p_initial + k * (lambda_d / (lambda_d - lambda_p)) * A_p_initial * (1 - exp(-(lambda_d - lambda_p) * t1))
term_k = k
term_lambda_ratio = lambda_d / (lambda_d - lambda_p)
term_decay_factor = lambda_d - lambda_p

print("Calculation Steps and Final Equation:")
print(f"1. From the measurements after separation (S(0)={S0:.1f}, S(14)={S14:.1f}), the relative Cherenkov efficiency 'k' was calculated to be: {k:.3f}")
print("\n2. The time before separation (t1) is found by solving the equation:")
print(f"Signal_unseparated(t1) = Signal_separated(14_days)")
print(f"A_p(t1) + k * A_d(t1) = {S14:.1f}")
print(f"{A_p_initial:.1f} + {term_k:.3f} * (({lambda_d:.4f} / ({lambda_d:.4f} - {lambda_p:.4f})) * {A_p_initial:.1f} * (1 - exp(-({lambda_d:.4f} - {lambda_p:.4f}) * t1))) = {S14:.1f}")
print(f"{A_p_initial:.1f} + {term_k:.3f} * ({term_lambda_ratio:.4f} * {A_p_initial:.1f} * (1 - exp(-{term_decay_factor:.4f} * t1))) = {S14:.1f}")

print(f"\nSolving this equation for t1 yields:")
print(f"Approximate time (t1): {t1:.2f} days")
print(f"This is approximately {t1*24:.1f} hours.")
print("\nFinal Answer:")
# Final answer is requested in days.
final_answer = round(t1, 2)
# To match the requested format, although printing it is better for the user
# print(f'<<<{final_answer}>>>')
print(f"{final_answer} days")