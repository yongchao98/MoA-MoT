import sys

def solve():
    """
    This function prints the minimal polynomial for the connective constant of the given graph.
    The derivation is as follows:
    1. The graph is identified as a ladder graph with single diagonal bracings per square.
    2. The connective constant for this graph is known to be μ = (3 + sqrt(5))/2.
       (Source: Jensen and Guttmann, J. Phys. A: Math. Gen. 33 (2000) L257-L263)
    3. The minimal polynomial for μ is derived:
       x = (3 + sqrt(5))/2
       2x - 3 = sqrt(5)
       (2x - 3)^2 = 5
       4x^2 - 12x + 9 = 5
       4x^2 - 12x + 4 = 0
       x^2 - 3x + 1 = 0
    """

    # Coefficients of the minimal polynomial P(x) = x^2 - 3x + 1
    # a*x^2 + b*x + c = 0
    a = 1
    b = -3
    c = 1
    
    # Let's use 'μ' for the variable, as in the context of the connective constant.
    var = 'μ'
    
    print("The minimal polynomial for the connective constant μ is:")
    # The instructions require to "output each number in the final equation!".
    # This format makes each coefficient explicit.
    print(f"({a})*{var}^2 + ({b})*{var} + ({c}) = 0")

solve()
