def calculate_profit_area():
    """
    This function calculates the firm's profit/loss using example values
    and demonstrates that it is represented by the area S - H.
    """
    # Let's assign some example values that satisfy the condition P < AVC < ATC
    # Let's assume the firm is producing at a quantity q1 = 100 units.
    q1 = 100
    P = 10   # The market price is $10.
    AVC = 12 # The average variable cost is $12.
    ATC = 15 # The average total cost is $15.

    # The condition P < AVC < ATC is met (10 < 12 < 15).

    # Step 1: Calculate the areas S, T, and H, which correspond to
    # Total Revenue, Total Variable Cost, and Total Cost, respectively.
    S = P * q1    # Area S = Total Revenue (TR)
    T = AVC * q1  # Area T = Total Variable Cost (TVC)
    H = ATC * q1  # Area H = Total Cost (TC)

    # Step 2: Calculate the firm's profit. Profit = Total Revenue - Total Cost.
    profit = S - H

    # Step 3: Print the breakdown of the calculation and the final equation.
    print("--- Economic Analysis with Example Values ---")
    print(f"Quantity (q1): {q1}")
    print(f"Price (P): {P}")
    print(f"Average Variable Cost (AVC): {AVC}")
    print(f"Average Total Cost (ATC): {ATC}\n")

    print("--- Calculating the Areas ---")
    print(f"Area S (Total Revenue) = P * q1 = {P} * {q1} = {S}")
    print(f"Area H (Total Cost)     = ATC * q1 = {ATC} * {q1} = {H}\n")

    print("--- Final Profit/Loss Calculation ---")
    print("The firm's profit is calculated as Total Revenue - Total Cost.")
    print("In terms of the defined areas, this corresponds to S - H.\n")

    print("Final Equation:")
    print(f"Profit = S - H")
    print(f"Profit = {S} - {H} = {profit}")
    
    # Since the profit is negative, the firm is incurring a loss of 500.
    # The area representing the profit is S - H.

calculate_profit_area()