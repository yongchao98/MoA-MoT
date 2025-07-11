def calculate_profit_area():
    """
    Calculates and explains the area representing a firm's profit or loss
    based on the given market conditions.
    """
    # 1. Define variables with example numbers that satisfy P < AVC < ATC.
    q1 = 50   # Output quantity
    P = 20    # Market Price, also equals MC
    AVC = 25  # Average Variable Cost
    ATC = 30  # Average Total Cost

    print("--- Problem Setup ---")
    print(f"A firm is producing at quantity q1 = {q1}.")
    print(f"The market condition is P < AVC < ATC.")
    print(f"Using example values: P={P}, AVC={AVC}, ATC={ATC}. The condition {P} < {AVC} < {ATC} is met.\n")

    # 2. Calculate the areas S and H.
    # Area S represents Total Revenue (TR = P * q1).
    S = P * q1
    # Area H represents Total Cost (TC = ATC * q1).
    H = ATC * q1

    print("--- Area Definitions ---")
    print(f"Area S represents Total Revenue (TR).")
    print(f"TR = P * q1 = {P} * {q1} = {S}")
    print(f"Area H represents Total Cost (TC).")
    print(f"TC = ATC * q1 = {ATC} * {q1} = {H}\n")

    # 3. Calculate the profit.
    # Profit (π) = Total Revenue (TR) - Total Cost (TC) = S - H.
    profit = S - H

    print("--- Profit Calculation ---")
    print("Profit (π) is defined as Total Revenue minus Total Cost (π = TR - TC).")
    print("In terms of the defined areas, this is π = S - H.")
    print("The area that represents the firm's profit is calculated as:")
    # The final output shows the equation with each number substituted.
    print(f"\nProfit = S - H = {S} - {H} = {profit}\n")
    print(f"Since the value is negative, the firm is experiencing a loss of {-profit}.")

calculate_profit_area()
