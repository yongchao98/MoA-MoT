import sympy

# Define symbols for the parameters and variables
theta, x = sympy.symbols('theta x')
t = sympy.Function('t')
p = sympy.Function('p')
F = sympy.Function('F')

# Define the cost function based on the problem description
# C(theta, x) = t(theta) + x*(1-theta)*p(x)*F(theta)
def cost_function(theta_val, x_val, t_func, p_func, F_func):
    """Calculates the expected cost."""
    # The term x*(1-theta) is the probability of refusal
    # p(x)*F(theta) is the expected penalty given refusal
    # t(theta) is the base tax liability
    return t_func(theta_val) + x_val * (1 - theta_val) * p_func(x_val) * F_func(theta_val)

# Define example decreasing functions for t, p, F.
# We assume non-negative probabilities and penalties.
# The exact values don't change the outcome, only their relationships (e.g., p(1)>0, F(0)>0).
t_vals = {0: 100, 1: 50}
p_vals = {0: 0.8, 1: 0.5} # p is decreasing, p(0) > p(1)
F_vals = {0: 1000, 1: 100} # F is decreasing, F(0) > F(1)

t_func = lambda val: t_vals[val]
p_func = lambda val: p_vals[val]
F_func = lambda val: F_vals[val]

# --- Step 1: Calculate costs for each type of firm and auditor ---
print("--- Calculating Expected Costs C(theta, x) ---")

# For malpractice firms (theta = 0)
c_0_0 = cost_function(0, 0, t_func, p_func, F_func)
c_0_1 = cost_function(0, 1, t_func, p_func, F_func)
print(f"Cost for malpractice firm (theta=0) with lenient auditor (x=0): C(0,0) = {t_vals[0]} + 0*(1-0)*p(0)*F(0) = {c_0_0}")
print(f"Cost for malpractice firm (theta=0) with strict auditor (x=1): C(0,1) = {t_vals[0]} + 1*(1-0)*{p_vals[1]}*{F_vals[0]} = {c_0_1}")

# For truthful firms (theta = 1)
c_1_0 = cost_function(1, 0, t_func, p_func, F_func)
c_1_1 = cost_function(1, 1, t_func, p_func, F_func)
print(f"Cost for truthful firm (theta=1) with lenient auditor (x=0): C(1,0) = {t_vals[1]} + 0*(1-1)*p(0)*F(1) = {c_1_0}")
print(f"Cost for truthful firm (theta=1) with strict auditor (x=1): C(1,1) = {t_vals[1]} + 1*(1-1)*{p_vals[1]}*{F_vals[1]} = {c_1_1}")
print("-" * 50)

# --- Step 2: Determine optimal auditor choice x*(theta) for each firm type ---
print("--- Determining Optimal Auditor Choice x*(theta) ---")
# For theta = 0
if c_0_0 < c_0_1:
    x_star_0 = 0
    print("For malpractice firm (theta=0), C(0,0) < C(0,1). The firm chooses the lenient auditor (x=0). So, x*(0) = 0.")
elif c_0_1 < c_0_0:
    x_star_0 = 1
    print("For malpractice firm (theta=0), C(0,1) < C(0,0). The firm chooses the strict auditor (x=1). So, x*(0) = 1.")
else:
    # Assuming tie-breaking towards leniency
    x_star_0 = 0
    print("For malpractice firm (theta=0), costs are equal. Assuming tie-break towards leniency, x*(0) = 0.")

# For theta = 1
if c_1_0 < c_1_1:
    x_star_1 = 0
elif c_1_1 < c_1_0:
    x_star_1 = 1
else:
    # A truthful firm is indifferent. We apply a reasonable tie-breaking rule.
    # Choosing the simpler option (lenient auditor) is a standard tie-breaking assumption.
    x_star_1 = 0
    print("For truthful firm (theta=1), C(1,0) = C(1,1). The firm is indifferent.")
    print("Assuming a tie-breaking rule to choose the more lenient auditor, x*(1) = 0.")
print("-" * 50)

# --- Step 3: Evaluate the condition and find the resulting set ---
print("--- Evaluating the Proposition ---")
# The proposition is: "companies keeping more truthful accounts choose more lenient auditors".
# This means that if theta_H > theta_L, then x*(theta_H) <= x*(theta_L).
# In our case, this means x*(1) <= x*(0).
proposition_holds = x_star_1 <= x_star_0

print(f"The proposition is equivalent to the condition x*(1) <= x*(0).")
print(f"With our calculated choices, this is: {x_star_1} <= {x_star_0}")
print(f"This condition is {proposition_holds}.")
print("\nSince the proposition is true for the model, it holds for any value of theta.")
print("The set of values of theta for which the proposition holds is the entire domain of theta.")

final_set = "{0, 1}"
print("\nThe final set of values for theta is:")
print(final_set)