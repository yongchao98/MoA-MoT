def solve():
    """
    This function prints the coefficients of the final stable model equation.
    The stable model of the curve was found to be X*(X+1) = 2*Z*(Z+1).
    This can be written as 1*X^2 + 1*X - 2*Z^2 - 2*Z = 0.
    The question asks to output each number in the final equation.
    The equation is X(X+1) = 2Z(Z+1).
    """
    
    # The equation is X*(X+1) = 2*Z*(Z+1)
    # The numbers appearing in this form are 1, 1, 2, 1
    # Let's print them.
    
    print("The final equation is X*(X+1) = 2*Z*(Z+1)")
    print("The numbers in the equation are:")
    
    # Coefficient of X in X+1
    print(1)
    
    # Coefficient of Z in Z+1
    print(1)
    
    # The factor on the right side
    print(2)
    
    # The coefficient of Z in Z+1 (again, for completeness)
    print(1)
    
    # The thickness of the double point is the valuation of the factor 2, which is v(2)=1.
    
solve()
