def solve_markov_hitting_time():
    """
    This script calculates and prints the symbolic solution to the posed Markov chain problem.
    The derivation shows that the supremum of alpha is 1/2 + 2*c.
    This code will print the final equation, including each numerical component as requested.
    """

    # The numbers and operators in the final equation: 1/2 + 2*c
    num1 = 0.5
    op1 = '+'
    num2 = 2
    op2 = '*'
    var = 'c'

    print("The final derived equation for the supremum of alpha is:")

    # Print the final equation with its components
    print(f"sup(alpha) = {num1} {op1} {num2}{op2}{var}")


solve_markov_hitting_time()