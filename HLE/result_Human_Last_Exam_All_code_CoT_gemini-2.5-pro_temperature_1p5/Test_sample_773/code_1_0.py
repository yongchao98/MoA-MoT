def solve_and_print_formula():
    """
    This function prints the derived formula for the total mass.
    The parameters n, q_v, and the function Z are treated as symbols as they
    are not specified in the problem statement.
    """

    # The final formula for the total mass is: (q_v / (q_v - 1)) * Product_{k=2 to n} Z(k)
    # The instruction is to "output each number in the final equation".
    # I will interpret this as explicitly stating the numerical coefficients and constants
    # that appear in the general formula.

    print("The total mass is given by a formula composed of two main parts.")
    print("\nPart 1: A factor dependent on q_v, the order of the residual field.")
    print("This factor is: q_v / (q_v - 1)")
    print("The numbers in this part of the equation are:")
    print("  - The coefficient of q_v in the numerator is: 1")
    print("  - The coefficient of q_v in the denominator is: 1")
    print("  - The constant term in the denominator is: -1")

    print("\nPart 2: A product of special values of the Dedekind zeta function Z.")
    print("This product is: Z(2) * Z(3) * ... * Z(n)")
    print("The numbers in this part of the equation are the arguments to the zeta function,")
    print("which are the integers from 2 to n, inclusive.")

    print("\n-------------------------------------------")
    print("The final formula for the total mass is:")
    print("Mass = (q_v / (q_v - 1)) * (Z(2) * Z(3) * ... * Z(n))")
    print("-------------------------------------------")

# Execute the function to print the solution.
solve_and_print_formula()