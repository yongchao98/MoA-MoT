import numpy as np

# Step 1: Define the coefficients of the cubic equation for the equilibrium quantity (Q_eq)
# The equation is 0.02*Q^3 - 0.09*Q^2 + 0*Q - 4 = 0
coeffs = [0.02, -0.09, 0, -4]

# Step 2: Solve for the roots of the equation
roots = np.roots(coeffs)

# Find the real root that is physically meaningful (positive)
# The producer's capacity is 10, so we look for a root in the range [0, 10]
Q_eq = 0
for r in roots:
    if np.isreal(r) and r > 0 and r <= 10:
        Q_eq = np.real(r)
        break

# Step 3: Calculate the excess demand using the derived formula
# Excess Demand = 45 * Q_eq^3
excess_demand = 45 * (Q_eq**3)

# Step 4: Print the results step-by-step
print("The producer maximizes profit where Marginal Revenue (MR) equals Marginal Cost (MC).")
print("With MC=0, the equation to solve is MR(Q) = 0.")
print(f"The equation for the equilibrium quantity (Q_eq) is: 0.02*Q_eq^3 - 0.09*Q_eq^2 - 4 = 0")
print(f"The profit-maximizing quantity supplied is Q_eq = {Q_eq:.4f}")
print("\nThe excess demand is calculated using the formula: Excess Demand = 45 * Q_eq^3")
print(f"Plugging in the value of Q_eq: Excess Demand = 45 * ({Q_eq:.4f})^3")
print(f"The value of the excess demand at the equilibrium price is: {excess_demand:.4f}")
