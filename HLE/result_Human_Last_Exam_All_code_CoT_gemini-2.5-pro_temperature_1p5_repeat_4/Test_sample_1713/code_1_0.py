def calculate_profit_area():
    """
    Calculates and explains the area representing a firm's profit or loss
    based on the given economic conditions.
    """
    # 1. Define variables with example values that fit the problem's condition:
    # P = MC < AVC < ATC
    q1 = 50   # Output level
    P = 20    # Market price (P), also equal to Marginal Cost (MC)
    AVC = 25  # Average Variable Cost
    ATC = 35  # Average Total Cost

    print("--- Firm Profit/Loss Calculation ---")
    print(f"Given conditions at output level q1 = {q1}:")
    print(f"Price (P) = {P}")
    print(f"Average Variable Cost (AVC) = {AVC}")
    print(f"Average Total Cost (ATC) = {ATC}")
    print("The condition P < AVC < ATC is satisfied.")
    print("-" * 35)

    # 2. Define profit, TR, and TC in terms of areas S and H.
    # Total Revenue (TR) is represented by area S.
    # TR = P * q1
    S = P * q1

    # Total Cost (TC) is represented by area H.
    # TC = ATC * q1
    H = ATC * q1

    # 3. The firm's profit or loss is the difference between TR and TC.
    profit_loss = S - H
    
    # 4. Explain the relationship and print the result.
    print("Step 1: The firm's Total Revenue (TR) is Price * quantity.")
    print(f"TR is represented by area S = P * q1 = {P} * {q1} = {S}\n")
    
    print("Step 2: The firm's Total Cost (TC) is Average Total Cost * quantity.")
    print(f"TC is represented by area H = ATC * q1 = {ATC} * {q1} = {H}\n")

    print("Step 3: The firm's profit or loss is TR - TC, which corresponds to S - H.")
    print("The area that represents the firm's profit or loss is therefore S - H.\n")
    
    print("--- Final Result ---")
    print(f"Profit/Loss = S - H = {S} - {H} = {profit_loss}")

    # Since P < ATC, the value is negative, indicating a loss.
    if profit_loss < 0:
        print(f"The firm is making a loss of {-profit_loss}.")
    else:
        print(f"The firm is making a profit of {profit_loss}.")

# Execute the function
calculate_profit_area()