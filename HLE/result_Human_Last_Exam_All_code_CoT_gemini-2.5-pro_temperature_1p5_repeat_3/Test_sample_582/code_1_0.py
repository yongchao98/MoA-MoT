def calculate_ex_dividend_price(q, E, d, g):
    """
    Calculates the per-share ex-dividend price in year 1 under the new policy.

    Args:
        q (float): Number of outstanding shares.
        E (float): Total market value of equity.
        d (float): Total dividends to be distributed in year 1 under the old policy.
        g (float): Annual dividend growth rate.
    """
    print("This script calculates the per-share ex-dividend price in year 1 under the new policy.")
    print("-----------------------------------------------------------------------------------")
    print("Inputs:")
    print(f"  - Initial outstanding shares (q): {q:,.0f}")
    print(f"  - Initial total market value of equity (E): ${E:,.2f}")
    print(f"  - Original year 1 total dividend (d): ${d:,.2f}")
    print(f"  - Annual dividend growth rate (g): {g:.2%}")
    print("-----------------------------------------------------------------------------------")

    # The formula for the per-share ex-dividend price in year 1 (p_1) is:
    # p_1 = (E * (1 + g) - d) / q

    # Calculate the numerator: E * (1 + g) - d
    numerator = E * (1 + g) - d
    
    # Calculate the per-share price
    p1 = numerator / q
    
    # Print the calculation steps
    print("Step-by-step calculation for the per-share ex-dividend price (p_1):")
    print(f"p_1 = (E * (1 + g) - d) / q")
    print(f"p_1 = ({E:,.2f} * (1 + {g}) - {d:,.2f}) / {q:,.0f}")
    print(f"p_1 = ({E * (1 + g):,.2f} - {d:,.2f}) / {q:,.0f}")
    print(f"p_1 = {numerator:,.2f} / {q:,.0f}")
    print("-----------------------------------------------------------------------------------")
    print(f"The final per-share ex-dividend price in year 1 is: ${p1:.2f}")

    # Returning the final answer in the required format
    return p1

# --- Example Data ---
# You can change these values to see how they affect the result.
q_shares = 2000000.0   # Number of outstanding shares
E_value = 100000000.0  # Total market value of equity
d_dividend = 4000000.0   # Total original dividends in year 1
g_growth = 0.06          # Annual dividend growth rate (6%)

# Run the calculation and store the final answer
final_answer = calculate_ex_dividend_price(q_shares, E_value, d_dividend, g_growth)

# The final answer is enclosed in <<< >>>
print(f"\n<<<{final_answer:.2f}>>>")
