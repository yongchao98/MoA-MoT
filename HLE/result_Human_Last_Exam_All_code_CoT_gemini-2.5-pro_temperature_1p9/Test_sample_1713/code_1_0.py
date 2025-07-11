import sys

def solve_firm_profit():
    """
    Calculates and explains the firm's profit or loss based on geometric areas.

    This function uses example values that satisfy the problem's condition:
    P = MC < AVC < ATC. It then calculates the areas representing Total
    Revenue (S) and Total Cost (H) to determine the firm's profit (S - H).
    """

    # Step 1: Define example values satisfying the condition P < AVC < ATC.
    # Let q_1 be the output level.
    q_1 = 100
    # Let P be the market price.
    P = 30
    # Let AVC be the Average Variable Cost at q_1.
    AVC = 40
    # Let ATC be the Average Total Cost at q_1.
    ATC = 55

    # Check if the values satisfy the condition P < AVC < ATC
    if not (P < AVC < ATC):
        print(f"Error: The chosen values P={P}, AVC={AVC}, ATC={ATC} do not satisfy the condition P < AVC < ATC.", file=sys.stderr)
        return

    # Step 2: Calculate the areas S and H.
    # S represents Total Revenue (TR = P * q_1).
    S = P * q_1

    # H represents Total Cost (TC = ATC * q_1).
    H = ATC * q_1

    # T is defined for context, representing Total Variable Cost (TVC = AVC * q_1)
    T = AVC * q_1
    
    # Step 3: Calculate the firm's profit or loss.
    # Profit = Total Revenue (TR) - Total Cost (TC)
    profit = S - H

    # Step 4: Print the explanation and the final result.
    print("Economic Analysis:")
    print("1. Profit for a firm is calculated as Total Revenue (TR) minus Total Cost (TC).")
    print("2. Total Revenue (TR) is Price × Quantity. In this problem, TR = P * q_1, which is the definition of area S.")
    print("3. Total Cost (TC) is Average Total Cost × Quantity. In this problem, TC = ATC * q_1, which is the definition of area H.")
    print("4. Therefore, the firm's profit is represented by the area S minus the area H.\n")

    print("Calculation with Example Values:")
    print(f"Quantity (q_1) = {q_1}")
    print(f"Price (P) = {P}")
    print(f"Average Total Cost (ATC) = {ATC}\n")
    
    print(f"Area S (Total Revenue) = P * q_1 = {P} * {q_1} = {S}")
    print(f"Area H (Total Cost)   = ATC * q_1 = {ATC} * {q_1} = {H}\n")

    print("The final formula for the area representing profit is S - H.")
    print("Using our values, the profit calculation is:")
    print(f"Profit = S - H = {S} - {H} = {profit}")
    print(f"\nSince the result is negative, the firm is making a loss. The magnitude of the loss is H - S = {H} - {S} = {-profit}.")

solve_firm_profit()