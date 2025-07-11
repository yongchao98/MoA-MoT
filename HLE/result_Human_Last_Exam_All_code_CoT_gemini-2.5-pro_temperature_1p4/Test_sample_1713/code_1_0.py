def calculate_profit_area():
    """
    Calculates and explains the firm's profit/loss based on the geometric areas S and H.
    """
    # Assigning example values that satisfy the condition P < AVC < ATC
    # Let's assume the firm produces 20 units.
    q1 = 20
    # Let the market price be $10.
    P = 10
    # Let the Average Variable Cost be $15.
    AVC = 15
    # Let the Average Total Cost be $18.
    ATC = 18

    # The problem conditions P < AVC < ATC are met: 10 < 15 < 18.

    print("Step 1: Define the economic concepts.")
    print("Profit (or Loss) = Total Revenue (TR) - Total Cost (TC)")
    print("-" * 30)

    print("Step 2: Relate the geometric areas to economic variables.")
    # Calculate Area S, which represents Total Revenue (TR).
    # S = P * q1
    S = P * q1
    print(f"Total Revenue (TR) is Price * Quantity.")
    print(f"The area S is defined by height P ({P}) and width q1 ({q1}).")
    print(f"Therefore, S = TR = {P} * {q1} = {S}")
    print("-" * 30)

    # Calculate Area H, which represents Total Cost (TC).
    # H = ATC * q1
    H = ATC * q1
    print(f"Total Cost (TC) is Average Total Cost * Quantity.")
    print(f"The area H is defined by height ATC ({ATC}) and width q1 ({q1}).")
    print(f"Therefore, H = TC = {ATC} * {q1} = {H}")
    print("-" * 30)

    print("Step 3: Calculate the profit or loss.")
    # Profit/Loss is TR - TC, which is S - H.
    profit_loss = S - H
    print("The firm's profit is given by the formula: Profit = S - H")
    print("\nUsing the example values:")
    # The final print statement shows the equation with the calculated numbers
    print(f"The firm's profit (loss) = {S} - {H} = {profit_loss}")
    print("\nSince the result is negative, the firm is experiencing a loss.")
    print("The area that represents the firm's loss is the difference between the area representing total cost (H) and the area representing total revenue (S).")

calculate_profit_area()