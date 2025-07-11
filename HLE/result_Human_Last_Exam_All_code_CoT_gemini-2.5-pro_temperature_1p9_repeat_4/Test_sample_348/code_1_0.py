import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

# Step 1: Define the supply and demand functions
def demand_p(q):
  """Defines the demand function P = 18 * e^(-arctan(Q))."""
  return 18 * np.exp(-np.arctan(q))

def supply_p(q):
  """Defines the supply function P = ln(Q^3 - 2)."""
  # The function is defined for Q^3 - 2 > 0
  with np.errstate(invalid='ignore'): # Suppress log of negative number warnings during solve
    return np.log(q**3 - 2)

# Step 2: Find the market equilibrium
def equilibrium_equation(q):
  """Equation to solve for equilibrium: Demand(Q) - Supply(Q) = 0."""
  return demand_p(q) - supply_p(q)

# An initial guess is needed for the numerical solver.
# Let's try Q=4: P_demand is ~4.78, P_supply is ~4.13. They are close.
initial_guess = 4.0
# fsolve finds the root of the equation, which is our equilibrium quantity Q_E
Q_E = fsolve(equilibrium_equation, initial_guess)[0]

# Calculate the equilibrium price P_E using the demand function
P_E = demand_p(Q_E)

print(f"Market Equilibrium found at:")
print(f"  - Equilibrium Quantity (Q_E): {Q_E:.4f}")
print(f"  - Equilibrium Price (P_E):    {P_E:.4f}\n")

# Step 3: Calculate Consumer Surplus (CS)
# CS = Integral from 0 to Q_E of (Demand(Q) - P_E) dQ
cs_integrand = lambda q: demand_p(q) - P_E
consumer_surplus, cs_error = quad(cs_integrand, 0, Q_E)

# Step 4: Calculate Producer Surplus (PS)
# The supply curve is defined for Q > 2^(1/3). This is the lower limit for integration.
q_min_supply = (2)**(1/3)
# PS = Integral from q_min_supply to Q_E of (P_E - Supply(Q)) dQ
ps_integrand = lambda q: P_E - supply_p(q)
producer_surplus, ps_error = quad(ps_integrand, q_min_supply, Q_E)

print(f"Calculating welfare components:")
print(f"  - Consumer Surplus (CS): {consumer_surplus:.4f}")
print(f"  - Producer Surplus (PS): {producer_surplus:.4f}\n")

# Step 5: Calculate Total Welfare
# Total Welfare (TW) = Consumer Surplus + Producer Surplus
total_welfare = consumer_surplus + producer_surplus

print("The final equation for total welfare is:")
print(f"Total Welfare = Consumer Surplus + Producer Surplus")
print(f"Total Welfare = {consumer_surplus:.4f} + {producer_surplus:.4f} = {total_welfare:.4f}")
