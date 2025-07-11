import math

def final_answer():
    """
    This function prints the derived exact formula for l(n) and the integer
    coefficients present in the formula, as requested.
    """
    
    # The derived formula is:
    # l(n) = (2 * (n^2 + 1) - 2 * (2*n - 1) * sqrt(n^2 - n + 1)) / n^2
    # We can write this in a general form:
    # l(n) = (c1*(n^c9 + c2) - c3*(c4*n + c5)*sqrt(c6*n^2 + c7*n + c8)) / n^c9
    
    c1 = 2
    c2 = 1
    c3 = 2
    c4 = 2
    c5 = -1
    c6 = 1
    c7 = -1
    c8 = 1
    c9 = 2
    
    # Printing the final formula as a string
    formula_str = f"({c1}*(n^{c9} + {c2}) - {c3}*({c4}*n - {-c5})*sqrt({c6}*n^2 - {-c7}*n + {c8})) / n^{c9}"
    
    print("The exact value of l(n) is given by the formula:")
    print(f"l(n) = {formula_str}")
    
    print("\nThe numbers in the final equation are:")
    print(c1)
    print(c2)
    print(c3)
    print(c4)
    print(c5)
    print(c6)
    print(c7)
    print(c8)
    print(c9)

if __name__ == '__main__':
    final_answer()
