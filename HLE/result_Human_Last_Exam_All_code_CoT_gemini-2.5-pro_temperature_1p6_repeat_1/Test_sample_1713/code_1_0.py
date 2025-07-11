def calculate_profit_area():
    """
    Calculates and explains the area representing a firm's profit or loss
    based on the given market conditions.
    """

    # Step 1: Define variables with example values that fit the condition
    # P = MC < AVC < ATC.
    q1 = 20  # Output quantity
    P = 10   # Market price (and Marginal Cost)
    AVC = 12 # Average Variable Cost
    ATC = 15 # Average Total Cost

    # Step 2: Calculate the areas S and H based on their definitions.
    # Area S = Total Revenue (TR) = Price * Quantity
    S = P * q1

    # Area H = Total Cost (TC) = Average Total Cost * Quantity
    H = ATC * q1

    # Step 3: Explain the relationship between the areas and economic concepts.
    print("--- Economic Analysis ---")
    print("1. The firm's Total Revenue (TR) is Price × Quantity.")
    print(f"   In terms of the given areas, TR is represented by Area S, which is P * q1 = {P} * {q1} = {S}.")
    print("\n2. The firm's Total Cost (TC) is Average Total Cost × Quantity.")
    print(f"   In terms of the given areas, TC is represented by Area H, which is ATC * q1 = {ATC} * {q1} = {H}.")
    print("\n3. A firm's profit or loss is calculated as Total Revenue - Total Cost (TR - TC).")
    print("   By substituting the areas, Profit/Loss = S - H.")

    # Step 4: Calculate the final profit/loss value and print the equation.
    profit_or_loss_value = S - H

    print("\n--- Final Calculation ---")
    print("Since P < ATC, the firm is operating at a loss.")
    print("The value representing the firm's loss is given by the equation S - H:")
    print(f"\nFinal Equation: {S} - {H} = {profit_or_loss_value}")


# Run the function to display the results.
calculate_profit_area()