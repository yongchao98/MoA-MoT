import numpy as np

def simulate_bachelier_step():
    """
    Demonstrates a single step of the Bachelier model, which is mathematically
    equivalent to Brownian motion.
    """
    # --- Parameters for the simulation ---
    # S0: Initial stock price
    S0 = 100.0
    # mu: Annual drift (the expected trend of the stock)
    mu = 5.0
    # sigma: Annual volatility (the amount of randomness)
    sigma = 20.0
    # T: Time horizon in years
    T = 1.0/252 # Simulate one trading day

    # The Bachelier model equation is: dS = mu*dt + sigma*dW_t
    # where dS is the change in price, dt is the change in time,
    # and dW_t is the random "Wiener process" or Brownian motion term.

    # --- Calculation ---
    # Calculate the drift component for one time step
    drift_component = mu * T

    # The random component is drawn from a normal distribution
    # Its standard deviation is sqrt(T)
    random_shock = np.random.normal(0, np.sqrt(T))
    stochastic_component = sigma * random_shock

    # Calculate the total change in price (dS)
    dS = drift_component + stochastic_component

    # Calculate the new stock price (S1)
    S1 = S0 + dS

    # --- Output the result ---
    print("Louis Bachelier modeled stock price changes using an equation similar to the one for Brownian Motion.")
    print("The basic formula is: New Price = Old Price + Drift + Random Shock\n")
    print("Let's simulate one step:\n")
    print(f"Old Price (S0): {S0:.4f}")
    print(f"Drift term (mu * dt): {drift_component:.4f}")
    print(f"Random Shock term (sigma * dW): {stochastic_component:.4f}\n")
    print("Final Equation with our simulated numbers:")
    print(f"{S1:.4f} = {S0:.4f} + {drift_component:.4f} + {stochastic_component:.4f}")

if __name__ == "__main__":
    simulate_bachelier_step()