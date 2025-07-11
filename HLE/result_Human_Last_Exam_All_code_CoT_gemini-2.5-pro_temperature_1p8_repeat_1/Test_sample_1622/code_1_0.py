import sys

def solve():
    """
    This function generates and prints the formula for P(n).
    L is a placeholder for ln(n).
    """

    # Coefficients for the term proportional to 1/n^2
    c22 = 3
    c21 = 2
    c20 = -2
    d2 = 24

    # Coefficients for the term proportional to 1/n^3
    c33 = 1
    c32 = 4
    c31 = -8
    c30 = 2
    d3 = 48
    
    # Building the string for the formula
    # Using replace for cleaner output, e.g., "+ -2" becomes "- 2"
    p2_numerator = f"({c22}*L**2 + {c21}*L {c20})".replace("+ -", "- ")
    p2_term = f"{p2_numerator}/({d2}*n**2)"

    p3_numerator = f"({c33}*L**3 + {c32}*L**2 {c31}*L + {c30})".replace("+ -", "- ")
    p3_term = f"{p3_numerator}/({d3}*n**3)"
    
    # We write P(n) in a way that makes its composition clear.
    final_formula = f"{p2_term} + {p3_term}"
    
    print(final_formula)

if __name__ == '__main__':
    solve()