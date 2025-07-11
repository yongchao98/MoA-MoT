def calculate_profit_or_loss():
    """
    Calculates and explains the firm's profit or loss based on geometric areas.

    The problem defines the market conditions for a firm:
    P = MC < AVC < ATC at an output level q1.

    It also defines areas based on these values:
    - S = P * q1  (Total Revenue)
    - T = AVC * q1 (Total Variable Cost)
    - H = ATC * q1 (Total Cost)

    The firm's profit (or loss) is calculated as Total Revenue - Total Cost.
    """

    # Assign example values that satisfy P < AVC < ATC
    q1 = 50  # Output level
    P = 10   # Market Price
    AVC = 12 # Average Variable Cost
    ATC = 15 # Average Total Cost

    # Calculate the areas S, T, and H
    S = P * q1
    T = AVC * q1
    H = ATC * q1

    # Calculate the profit or loss
    # Profit = Total Revenue (S) - Total Cost (H)
    profit = S - H

    print("Step 1: Define the economic formulas.")
    print("Profit = Total Revenue (TR) - Total Cost (TC)")
    print("-" * 30)

    print("Step 2: Relate the formulas to the given areas S and H.")
    print(f"Total Revenue (TR) = Price * Quantity = {P} * {q1} = {S}. This is area S.")
    print(f"Total Cost (TC) = Average Total Cost * Quantity = {ATC} * {q1} = {H}. This is area H.")
    print("-" * 30)

    print("Step 3: Calculate the profit using the areas S and H.")
    print("The formula for profit is: Profit = S - H")
    print("\nUsing the example values:")
    # The final equation with numbers, as requested.
    print(f"Profit = {S} - {H} = {profit}")
    print("\nSince the result is negative, the firm is making a loss.")
    print(f"The area representing the magnitude of the loss is H - S = {H} - {S} = {-profit}.")

calculate_profit_or_loss()