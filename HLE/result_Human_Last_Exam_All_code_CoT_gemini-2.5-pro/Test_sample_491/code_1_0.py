import sympy as sp

# Plan:
# 1. Define the symbolic variable 'd'.
# 2. Determine Player 1's payoff in the standard Subgame Perfect Equilibrium (SPE).
# 3. Analyze Player 2's unique strategy to determine Player 1's optimal payoff against them.
# 4. Calculate the difference between Player 1's optimal payoff and the SPE payoff.
# 5. Express this difference as a simplified rational function p(d)/q(d).
# 6. Extract the integer coefficients from both the numerator p(d) and the denominator q(d).
# 7. Sum the absolute values of these coefficients to get the final answer.

# Step 1: Define the symbolic variable
d = sp.Symbol('d')

# Step 2: Player 1's payoff in the standard SPE
V1_spe = 1 / (1 + d)

# Step 3: Player 1's optimal payoff against the specific Player 2
# Let V_sim be the expected payoff for a player in a simulation where the pie size is 1.
# The simulation rules are symmetric.
# The equation for V_sim is: V_sim = (prob_accept)*E[payoff_accept] + (prob_reject)*E[payoff_reject]
# V_sim = 0.5 * 0.5 + 0.5 * d * V_sim
V_sim_symbol = sp.Symbol('V_sim')
eq_V_sim = sp.Eq(V_sim_symbol, 0.25 + 0.5 * d * V_sim_symbol)
v_sim_sol = sp.solve(eq_V_sim, V_sim_symbol)[0]

# Player 2's rejection threshold is the expected payoff from an R-simulation.
# This simulation starts after rejection, so pie is d and it's P2's turn to propose.
p2_rejection_payoff = d * v_sim_sol

# Player 1's optimal payoff is by offering P2 an amount just above this threshold.
V1_optimal = 1 - p2_rejection_payoff

# Step 4: Calculate the difference
difference = V1_optimal - V1_spe

# Step 5: Simplify to a rational function p(d)/q(d)
diff_fraction = sp.cancel(difference)
p_d, q_d = sp.fraction(diff_fraction)

# Step 6: Extract integer coefficients
# For a canonical form, we make the leading coefficient of the denominator positive.
if q_d.as_poly(d).LC() < 0:
    p_d = -p_d
    q_d = -q_d

p_poly = sp.Poly(p_d, d)
q_poly = sp.Poly(q_d, d)

p_coeffs = [int(c) for c in p_poly.all_coeffs()]
q_coeffs = [int(c) for c in q_poly.all_coeffs()]

# Step 7: Sum the absolute values of the coefficients
sum_abs_p = sum(abs(c) for c in p_coeffs)
sum_abs_q = sum(abs(c) for c in q_coeffs)
total_sum = sum_abs_p + sum_abs_q

# Output the numbers in the final equation and the result
p_str = " + ".join(f"({c})*d**{p_poly.degree()-i}" for i, c in enumerate(p_coeffs)).replace("+ (-", "- (")
q_str = " + ".join(f"({c})*d**{q_poly.degree()-i}" for i, c in enumerate(q_coeffs)).replace("+ (-", "- (")

print(f"The difference is p(d)/q(d).")
print(f"Numerator p(d) = {p_str}")
print(f"Denominator q(d) = {q_str}")
print(f"")
print(f"The integer coefficients of the numerator p(d) are: {p_coeffs}")
print(f"The integer coefficients of the denominator q(d) are: {q_coeffs}")
print(f"")
print(f"The sum of the absolute values of coefficients in p(d) is: {sum_abs_p}")
print(f"The sum of the absolute values of coefficients in q(d) is: {sum_abs_q}")
print(f"The total sum of the absolute values of all coefficients is: {total_sum}")
<<<14>>>