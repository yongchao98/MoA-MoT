import math

def calculate_feynman_diagrams(n):
    """
    Calculates a(n), the number of non-vanishing Feynman diagrams of order 2n.
    The formula is a(n) = (2n-1)!!, where !! is the double factorial.
    """
    # The value to calculate the double factorial for is k = 2n - 1
    k = 2 * n - 1

    print(f"The number of non-vanishing Feynman diagrams of order 2n is given by the formula a(n) = (2n-1)!!")
    print(f"We want to find a(n) for n = {n}.")
    print(f"First, we calculate the term inside the double factorial: 2*n - 1 = 2*{n} - 1 = {k}")
    print(f"So, we need to calculate a({n}) = {k}!!")

    # Calculate the double factorial (2n-1)!!
    # The double factorial k!! for an odd k is the product of all odd integers from 1 to k.
    result = 1
    equation_terms = []
    for i in range(k, 0, -2):
        result *= i
        equation_terms.append(str(i))
    
    equation_str = " * ".join(equation_terms)
    
    print(f"The double factorial {k}!! is the product of all odd integers from 1 up to {k}.")
    print(f"So, the final equation is a({n}) = {equation_str} = {result}")

# The user is asking for a(3).
n_value = 3
calculate_feynman_diagrams(n_value)