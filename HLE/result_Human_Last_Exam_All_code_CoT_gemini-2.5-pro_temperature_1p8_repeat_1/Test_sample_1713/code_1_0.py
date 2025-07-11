def calculate_profit_loss_area():
    """
    This function demonstrates the calculation of a firm's loss in a perfectly competitive market
    under the condition P = MC < AVC < ATC.
    """
    # 1. Assign example values that satisfy the condition P < AVC < ATC.
    # We choose an arbitrary quantity q_1 and corresponding price and costs.
    q_1 = 20  # Output quantity
    P = 50    # Market Price
    AVC_q1 = 60 # Average Variable Cost at q_1
    ATC_q1 = 75 # Average Total Cost at q_1

    print(f"Let's analyze the firm's situation at output quantity q_1 = {q_1}.")
    print(f"The given market price and cost structure is:")
    print(f"Price (P) = ${P}")
    print(f"Average Variable Cost (AVC) = ${AVC_q1}")
    print(f"Average Total Cost (ATC) = ${ATC_q1}")
    print(f"This satisfies the problem's condition P < AVC < ATC ({P} < {AVC_q1} < {ATC_q1}).\n")

    # 2. Calculate the areas S (Total Revenue) and H (Total Cost).
    # Area S is a rectangle with height P and width q_1.
    S = P * q_1
    # Area H is a rectangle with height ATC and width q_1.
    H = ATC_q1 * q_1

    print("Step 1: Calculate Total Revenue (TR) and Total Cost (TC).")
    print(f"Total Revenue (TR) is represented by Area S = P * q_1 = {P} * {q_1} = {S}.")
    print(f"Total Cost (TC) is represented by Area H = ATC * q_1 = {ATC_q1} * {q_1} = {H}.\n")

    # 3. Calculate the profit or loss using the formula: Profit = TR - TC = S - H.
    profit_or_loss = S - H
    print("Step 2: Calculate the firm's profit or loss.")
    print(f"Profit = TR - TC = S - H")
    print(f"Profit = {S} - {H} = {profit_or_loss}")

    # 4. Determine the area representing the loss.
    # Since profit is negative, the firm has a loss. The area of the loss is H - S.
    if profit_or_loss < 0:
        loss_amount = -profit_or_loss
        loss_area = H - S
        print(f"The negative result indicates a loss of ${loss_amount}.\n")
        print("Step 3: Determine the area that represents this loss.")
        print("The area of the loss is the difference between the Total Cost area (H) and the Total Revenue area (S).")
        print(f"Area of Loss = H - S = {H} - {S} = {loss_area}")
        print("\nConclusion: The area representing the firm's loss is H - S.")
    else:
        # This case won't happen under the problem's conditions.
        print("The firm is making a profit, represented by the area S - H.")

calculate_profit_loss_area()