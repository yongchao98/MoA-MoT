import math

def solve():
    """
    Calculates the highest possible order for the inertial quotient E.
    This is equivalent to calculating the order of the general linear group GL_4(2).
    """
    n = 4
    q = 2

    # The order of GL(n, q) is the product of (q^n - q^i) for i from 0 to n-1.
    terms = []
    for i in range(n):
        term = q**n - q**i
        terms.append(term)

    # Calculate the final result
    result = 1
    for term in terms:
        result *= term

    # Format the output string
    term_strings = [str(t) for t in terms]
    equation_str = " * ".join(term_strings) + " = " + str(result)
    
    print("The highest order that E can have is the order of Aut(D), which is |GL_4(2)|.")
    print("The calculation is the product of terms (2^4 - 2^i) for i=0, 1, 2, 3.")
    print("The final equation showing the calculation is:")
    print(equation_str)

solve()