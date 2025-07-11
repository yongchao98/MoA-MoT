import sympy

def display_formula():
    """
    This function explains and prints the formula for the second voltage plateau
    of a graphite anode during lithium intercalation.
    """
    
    # Define symbolic representations for the components of the formula
    V_plateau_2 = sympy.Symbol("V_plateau_2")
    mu_2 = sympy.Symbol("μ_2")
    mu_3 = sympy.Symbol("μ_3")
    e = sympy.Symbol("e")

    # The formula derived from the physical model
    formula = sympy.Eq(V_plateau_2, (mu_2 - mu_3) / e)
    
    # The problem asks for the formula expression itself.
    expression = formula.rhs

    print("The simple formula that best approximates the second voltage plateau is:")
    
    # Print the full equation in a readable format
    # The prompt requires showing the numbers in the equation, which are the indices '2' and '3'.
    print(f"V = ({mu_2} - {mu_3}) / {e}")
    
if __name__ == "__main__":
    display_formula()