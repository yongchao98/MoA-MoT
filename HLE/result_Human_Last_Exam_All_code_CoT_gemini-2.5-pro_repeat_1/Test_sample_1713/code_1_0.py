def calculate_loss_area():
    """
    This function demonstrates the calculation of the firm's loss area
    based on the principles of microeconomics.
    """
    # Step 1: Define variables with example values that satisfy the condition
    # P = MC < AVC < ATC. Let's assume the firm produces 20 units.
    q_1 = 20  # Output level
    P = 10    # Market Price (and Marginal Cost)
    AVC = 12  # Average Variable Cost
    ATC = 15  # Average Total Cost

    print(f"Given the economic conditions at output level q_1 = {q_1}:")
    print(f"Price (P) = {P}")
    print(f"Average Total Cost (ATC) = {ATC}")
    print(f"The condition P < ATC ({P} < {15}) is met, so the firm is making a loss.")
    print("-" * 30)

    # Step 2: Calculate the areas S (Total Revenue) and H (Total Cost).
    # S represents Total Revenue (TR = P * q_1)
    S = P * q_1
    # H represents Total Cost (TC = ATC * q_1)
    H = ATC * q_1

    print("Step 1: Calculate Total Revenue (Area S)")
    print(f"Total Revenue (S) = P * q_1 = {P} * {q_1} = {S}")
    print("\nStep 2: Calculate Total Cost (Area H)")
    print(f"Total Cost (H) = ATC * q_1 = {ATC} * {q_1} = {H}")
    print("-" * 30)

    # Step 3: The area representing the firm's loss is the difference
    # between Total Cost (H) and Total Revenue (S).
    loss_area = H - S

    print("Step 3: Calculate the area representing the firm's loss.")
    print("The area of the loss is Total Cost - Total Revenue.")
    print("\nFinal Equation:")
    print(f"Loss Area = H - S = {H} - {S} = {loss_area}")

calculate_loss_area()