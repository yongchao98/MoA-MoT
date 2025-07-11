import numpy as np

# Step 1 & 2: Derive demand, profit, and marginal revenue.
# Individual demand: q_i = 400 - 100P + Q/100 + 3Q^2 - Q^3/20
# Market demand Q = 100 * q_i => q_i = Q/100
# Substitute q_i into its own equation:
# Q/100 = 400 - 100P + (Q/100) + 3Q^2 - Q^3/20 (mistake here in problem description, this doesn't work, assuming Q = 100 * q_i holds for the entire expression)
# Correct derivation: Q = 100 * (400 - 100P + Q/100 + 3Q^2 - Q^3/20)
# Q = 40000 - 10000P + Q + 300Q^2 - 5Q^3
# 10000P = 40000 + 300Q^2 - 5Q^3
# Inverse demand P(Q) = 4 + 0.03*Q^2 - 0.0005*Q^3
# Profit pi(Q) = P(Q)*Q = 4Q + 0.03*Q^3 - 0.0005*Q^4 for Q <= 10
# Marginal Revenue MR(Q) = d(pi)/dQ = 4 + 0.09*Q^2 - 0.002*Q^3

# Step 3: Determine the Optimal Quantity Supplied (Q_s)
# Producer's MC=0 for Q<=10. Profit is maximized when MR=MC=0 or at the boundary (Q=10).
# Let's check MR at the boundary Q=10.
q_s = 10
mr_at_10 = 4 + 0.09 * (q_s**2) - 0.002 * (q_s**3)
# MR(Q) = 4 + 0.09*Q^2 - 0.002*Q^3.
# The derivative d(MR)/dQ = 0.18*Q - 0.006*Q^2 = 0.006Q(30-Q), which is positive for 0<Q<30.
# Since MR(0)=4 and MR is increasing on [0, 10], MR is always positive on [0, 10].
# This means profit is always increasing for Q in [0, 10].
# Therefore, the producer chooses the maximum possible quantity.
print(f"The marginal revenue at Q=10 is {mr_at_10}, which is positive.")
print("Since MR > MC (0) for the entire production range [0, 10], the producer will supply the maximum possible quantity.")
q_s = 10
print(f"Quantity Supplied (Q_s) = {q_s}\n")

# Step 4: Determine the Equilibrium Price (P_e)
# Price is set to sell Q_s=10 units. P_e = P(10)
p_e = 4 + 0.03 * (q_s**2) - 0.0005 * (q_s**3)
print(f"To sell {q_s} units, the producer sets the price P_e = {p_e}\n")

# Step 5 & 6: Find the stable quantity demanded (Q_d) at P_e
# At P_e=6.5, the market demand Q_d must satisfy:
# 10000 * P_e = 40000 + 300*Q_d^2 - 5*Q_d^3
# 10000 * 6.5 = 40000 + 300*Q_d^2 - 5*Q_d^3
# 65000 = 40000 + 300*Q_d^2 - 5*Q_d^3
# 5*Q_d^3 - 300*Q_d^2 + 25000 = 0
# Q_d^3 - 60*Q_d^2 + 5000 = 0
# We need to find the roots of this cubic equation.
coeffs_demand = [1, -60, 0, 5000]
roots = np.roots(coeffs_demand)

print(f"At price P_e = {p_e}, the possible quantities demanded are the roots of the equation Q^3 - 60Q^2 + 5000 = 0.")
print("The roots are:", [f"{r:.4f}" for r in roots])

# The roots are Q_d=10, Q_d=58.54, and Q_d=-8.54.
# We discard the negative root.
# Due to the network externalities, the equilibrium at Q_d=10 is unstable, while the one at Q_d=58.54 is stable.
# The market will converge to the stable equilibrium.
q_d_stable = max(r for r in roots if r > 0)
print(f"The stable equilibrium for quantity demanded is Q_d = {q_d_stable:.4f}\n")

# Step 7: Calculate Excess Demand
excess_demand = q_d_stable - q_s
print("The final calculation for excess demand is:")
print(f"Excess Demand = (Stable Quantity Demanded) - (Quantity Supplied)")
print(f"Excess Demand = {q_d_stable:.4f} - {q_s} = {excess_demand:.4f}")