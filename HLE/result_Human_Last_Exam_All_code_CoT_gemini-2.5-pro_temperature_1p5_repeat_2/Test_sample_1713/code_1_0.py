def calculate_profit_expression():
    """
    This function analyzes the firm's profit/loss based on the provided areas S, T, and H.
    It uses an example to illustrate the calculation.
    """
    
    # 1. Assign example values that satisfy the condition P < AVC < ATC.
    q1 = 100  # Output level
    P = 10    # Market Price, also equals Marginal Cost (MC)
    AVC = 12  # Average Variable Cost
    ATC = 15  # Average Total Cost

    # 2. Calculate the areas S and H based on their definitions.
    # S is the area representing Total Revenue (TR).
    S = P * q1
    
    # H is the area representing Total Cost (TC).
    H = ATC * q1

    # 3. Calculate the firm's profit.
    # Profit is defined as Total Revenue (TR) minus Total Cost (TC).
    # In terms of the given areas, this is S - H.
    profit_or_loss = S - H

    # 4. Print the explanation and the final calculation.
    print("Step 1: Define Profit")
    print("A firm's profit is calculated as Total Revenue (TR) minus Total Cost (TC).")
    print("Profit = TR - TC\n")
    
    print("Step 2: Relate Areas to Economic Concepts")
    print(f"Area S is a rectangle with height P ({P}) and width q1 ({q1}).")
    print(f"Its area is P * q1, which is the definition of Total Revenue (TR). So, S = TR.")
    print(f"Area H is a rectangle with height ATC ({ATC}) and width q1 ({q1}).")
    print(f"Its area is ATC * q1, which is the definition of Total Cost (TC). So, H = TC.\n")

    print("Step 3: Derive the Final Expression")
    print("Substituting S for TR and H for TC into the profit formula, we get:")
    print("Profit = S - H\n")
    
    print("Step 4: Calculate Using Example Values")
    print(f"Using the example values (q1={q1}, P={P}, ATC={ATC}):")
    print(f"Area S (TR) = {P} * {q1} = {S}")
    print(f"Area H (TC) = {ATC} * {q1} = {H}\n")

    print("The final equation for the firm's profit/loss is:")
    # The prompt requires outputting each number in the final equation.
    print(f"Profit/Loss = S - H = {S} - {H} = {profit_or_loss}")
    
    print("\nSince the result is negative, the firm is experiencing a loss.")
    print("The expression representing the firm's profit (or loss) is S - H.")

# Execute the function to display the result.
calculate_profit_expression()