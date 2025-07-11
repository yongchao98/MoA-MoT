import sympy
from sympy import S, Sum, oo, factorial, exp, Symbol

def solve_expected_value():
    """
    This function calculates the expected value of T by symbolically summing an infinite series.
    The expected value E[T] is given by the series:
    E[T] = sum_{i=1 to infinity} [ i / (3^i * (i+1)!) ]
    """
    
    # Define the index of summation
    i = Symbol('i', integer=True, positive=True)

    # Define the general term of the infinite series
    term = i / (3**i * factorial(i+1))
    
    # Create the symbolic representation of the infinite sum
    infinite_sum = Sum(term, (i, 1, oo))

    # Use sympy's doit() method to evaluate the sum
    result = infinite_sum.doit()

    # The result is of the form A + B * exp(C). We will display these components.
    # From mathematical derivation, we know A=3, B=-2, C=1/3.
    A = 3
    B = -2
    C = S(1)/3 # Using sympy's rational to maintain precision and nice printing

    print("The expected value E[T] is calculated by summing an infinite series.")
    print("The result is of the form: A + B * exp(C)")
    print("\nThe components of this expression are:")
    print(f"A = {A}")
    print(f"B = {B}")
    print(f"C = {C}")

    print("\nThe final equation for the expected value is:")
    print(f"E[T] = {A} + ({B}) * exp({C})")

    print("\nFor verification, the symbolic result computed by sympy is:")
    print(result)
    
    print("\nAnd its numerical approximation is:")
    print(result.evalf())

if __name__ == '__main__':
    solve_expected_value()