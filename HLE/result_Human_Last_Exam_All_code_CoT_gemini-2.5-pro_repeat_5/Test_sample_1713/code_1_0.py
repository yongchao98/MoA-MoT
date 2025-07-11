def calculate_profit_area():
    """
    This function determines and prints the expression for the firm's profit or loss
    based on the areas S and H, which represent Total Revenue and Total Cost, respectively.
    """

    # The firm's profit (or loss) is Total Revenue (TR) minus Total Cost (TC).
    # From the problem description:
    # Total Revenue (TR) = Price * Quantity = P * q1 = Area S
    # Total Cost (TC) = Average Total Cost * Quantity = ATC * q1 = Area H
    
    # Therefore, Profit = TR - TC = S - H.

    # We define the parts of the final equation.
    total_revenue_area = "S"
    operator = "-"
    total_cost_area = "H"

    # Print the final equation representing the firm's profit or loss.
    print("The area representing the firm's profit or loss is expressed as:")
    print(f"{total_revenue_area} {operator} {total_cost_area}")

calculate_profit_area()