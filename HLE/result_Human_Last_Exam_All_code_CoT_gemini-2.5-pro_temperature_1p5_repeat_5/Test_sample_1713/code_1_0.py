def calculate_loss_area():
    """
    Calculates and explains the area representing a firm's loss based on given economic conditions.

    The problem states: P = MC < AVC < ATC.
    We need to find the area representing the firm's profit or loss.

    1. Profit = Total Revenue (TR) - Total Cost (TC)
    2. TR = Price * Quantity = P * q1. This corresponds to area S.
    3. TC = Average Total Cost * Quantity = ATC * q1. This corresponds to area H.
    4. Therefore, Profit = S - H.
    5. Since P < ATC, it follows that S < H, meaning the firm is making a loss.
    6. The area representing the magnitude of this loss is TC - TR, which is H - S.
    """

    # Assign illustrative values that satisfy the condition P < AVC < ATC
    q1 = 50   # Quantity
    P = 10    # Price
    AVC = 12  # Average Variable Cost
    ATC = 15  # Average Total Cost

    print(f"Let's use some example values that fit the criteria (P < ATC):")
    print(f"Quantity (q1) = {q1}")
    print(f"Price (P) = ${P}")
    print(f"Average Total Cost (ATC) = ${ATC}\n")

    # Calculate the areas S and H as defined in the problem
    S = P * q1      # Represents Total Revenue (TR)
    H = ATC * q1    # Represents Total Cost (TC)

    print(f"First, we calculate the area S (Total Revenue):")
    print(f"S = P * q1 = {P} * {q1} = {S}\n")

    print(f"Next, we calculate the area H (Total Cost):")
    print(f"H = ATC * q1 = {ATC} * {q1} = {H}\n")

    # The profit/loss is TR - TC, which is S - H.
    # The area representing the loss is TC - TR, which is H - S.
    loss_area = H - S

    print("The profit is S - H. Since S < H, the firm has a loss.")
    print("The area that represents the magnitude of this loss is H - S.")
    print("\nFinal Equation:")
    print(f"Loss Area = H - S = {H} - {S} = {loss_area}")


calculate_loss_area()