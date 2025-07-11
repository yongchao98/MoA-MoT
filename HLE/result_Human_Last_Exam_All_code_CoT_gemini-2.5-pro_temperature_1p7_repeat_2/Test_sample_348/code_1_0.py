import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

# Step 1: Define supply and demand functions
def supply(Q_s):
    """
    Supply function P = ln(Q_s^3 - 2).
    Returns a very large number if Q_s is outside the domain to guide the solver.
    """
    if np.any(Q_s**3 <= 2):
        return np.inf
    return np.log(Q_s**3 - 2)

def demand(Q_d):
    """
    Demand function P = 18 * e^(-arctan(Q_d)).
    """
    return 18 * np.exp(-np.arctan(Q_d))

# Step 2: Find the market equilibrium point (Qe, Pe)
def equilibrium_equation(Q):
    """
    The equation to solve: Demand(Q) - Supply(Q) = 0.
    """
    return demand(Q) - supply(Q)

# Numerically solve for the equilibrium quantity (Qe)
# An initial guess of 4.0 is based on preliminary analysis of the functions
initial_guess_Q = 4.0
Q_e = fsolve(equilibrium_equation, initial_guess_Q)[0]

# Calculate the equilibrium price (Pe) using the demand function
P_e = demand(Q_e)

print("--- Market Equilibrium ---")
print(f"The equilibrium equation to solve is: 18 * exp(-arctan(Q)) = ln(Q^3 - 2)")
print(f"Equilibrium Quantity (Qe): {Q_e:.4f}")
print(f"Equilibrium Price (Pe): {P_e:.4f}")
print("-" * 26)

# Step 3: Calculate Consumer Surplus (CS)
# CS is the integral of (Demand(Q) - Pe) from 0 to Qe
cs_integrand = lambda Q: demand(Q) - P_e
consumer_surplus, _ = quad(cs_integrand, 0, Q_e)

print("\n--- Consumer and Producer Surplus ---")
print(f"Consumer Surplus (CS) is the integral of (Demand(Q) - Pe) from Q=0 to Q={Q_e:.4f}.")
print(f"CS = ∫[0 to {Q_e:.4f}] (18 * e^(-arctan(Q)) - {P_e:.4f}) dQ = {consumer_surplus:.4f}")


# Step 4: Calculate Producer Surplus (PS)
# The supply curve starts at Q = 2^(1/3)
q_min_supply = 2**(1/3)
# PS is the integral of (Pe - Supply(Q)) from q_min_supply to Qe
ps_integrand = lambda Q: P_e - supply(Q)
producer_surplus, _ = quad(ps_integrand, q_min_supply, Q_e)

print(f"\nProducer Surplus (PS) is the integral of (Pe - Supply(Q)) from Q={q_min_supply:.4f} to Q={Q_e:.4f}.")
print(f"PS = ∫[{q_min_supply:.4f} to {Q_e:.4f}] ({P_e:.4f} - ln(Q^3 - 2)) dQ = {producer_surplus:.4f}")
print("-" * 37)

# Step 5: Calculate Total Welfare
total_welfare = consumer_surplus + producer_surplus

print("\n--- Total Welfare ---")
print(f"Total Welfare = Consumer Surplus + Producer Surplus")
print(f"Total Welfare = {consumer_surplus:.4f} + {producer_surplus:.4f}")
print(f"Final Total Welfare: {total_welfare:.4f}")

<<<11.5204>>>