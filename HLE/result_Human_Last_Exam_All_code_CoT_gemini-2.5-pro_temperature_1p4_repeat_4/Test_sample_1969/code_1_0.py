import sympy

def solve():
    """
    This function derives and prints the formula for c_k, the Shapley value for player p_k.
    """
    n, k = sympy.symbols('n k')

    # The simplified formula for the Shapley value c_k is:
    # c_k = k*S1**3 - k**2*S1**2 + k*S1*S2
    # where S1 is the sum of the first n integers and S2 is the sum of the first n squares.

    S1 = n * (n + 1) / 2
    S2 = n * (n + 1) * (2 * n + 1) / 6

    # Substitute S1 and S2 into the formula for c_k
    c_k_formula = k * S1**3 - k**2 * S1**2 + k * S1 * S2

    # Simplify the expression
    c_k_simplified = sympy.factor(c_k_formula)
    
    # The prompt requests the formula be printed clearly, showing all the numbers.
    # The simplified formula is c_k = (k * n**2 * (n + 1)**2 * ((n + 2)*(3*n + 1) - 6*k)) / 24.
    # We will format the print statement to reflect this structure.

    # Using pretty print for a more readable mathematical formula output
    # but the request asks to print numbers, so a f-string is better
    
    print("The formula for the exact amount of money c_k that player p_k gets is:")
    
    # We build the string manually to match the required formatting style and include all numbers.
    # The structure is c_k = k/24 * n^2 * (n+1)^2 * [ ... ]
    
    term1 = "k"
    term2 = "n**2"
    term3 = "(n + 1)**2"
    term4 = "((n + 2)*(3*n + 1) - 6*k)"
    denominator = "24"
    
    print(f"c_k = ({term1} * {term2} * {term3} * {term4}) / {denominator}")
    

solve()