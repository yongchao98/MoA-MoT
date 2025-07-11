import sympy

def solve_expected_value():
    """
    Calculates the expected value of T using symbolic mathematics.
    The derivation leads to an infinite series, which is computed here using sympy.
    """
    # Define the symbolic variable for the summation index
    i = sympy.Symbol('i', integer=True, positive=True)

    # As derived in the thinking steps, the expected value E[T] is given by the sum:
    # E[T] = Sum_{i=1 to infinity} 1 / (3^i * (i+1) * (i-1)!)
    # To compute this sum, we can rewrite the general term.
    # 1 / (3^i * (i+1) * (i-1)!) is equivalent to i / (3**i * (i+1)!).
    term = i / (sympy.factorial(i + 1) * 3**i)

    # We use sympy's summation function to find the exact value of the infinite series.
    expected_value = sympy.Sum(term, (i, 1, sympy.oo)).doit()

    # The computed result is 3 - 2*exp(1/3). We will print this result.
    print(f"The exact expected value of T is: {expected_value}")

    # The problem requires printing each number in the final equation.
    # The final equation is E[T] = 3 - 2 * exp(1/3).
    # The numbers involved in this expression are 3, 2, 1, and 3.
    print("The numbers appearing in the final expression 3 - 2 * exp(1/3) are:")
    print(3)
    print(2)
    print(1)
    print(3)

solve_expected_value()