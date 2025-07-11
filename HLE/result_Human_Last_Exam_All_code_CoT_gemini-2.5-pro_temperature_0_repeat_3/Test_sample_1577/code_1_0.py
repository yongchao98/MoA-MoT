def solve_toric_code_gsd():
    """
    Calculates and prints the formula for the ground space degeneracy (GSD)
    of the toric code on a sphere with n smooth and m rough holes.

    The GSD is 2^k, where k is the number of logical qubits.
    The formula for k on a surface with genus g, n smooth boundaries,
    and m rough boundaries is:
        k = 2*g + n + m - c
    where c is the number of boundary types present.

    For a sphere, g=0. In the general case with both smooth (n>0) and
    rough (m>0) holes, there are c=2 types of boundaries.
    So, k = 2*0 + n + m - 2 = n + m - 2.

    The GSD is therefore 2^(n + m - 2).
    """
    # The formula is derived from topological quantum field theory principles.
    # The variables n and m represent the number of smooth and rough holes, respectively.
    base = 2
    exponent_term_n = "n"
    exponent_term_m = "m"
    exponent_constant = -2

    # We will print the formula as a string.
    # The problem asks to output each number in the final equation.
    # The numbers are the base '2' and the constant '-2' in the exponent.
    print(f"The ground space degeneracy (GSD) is given by the formula:")
    print(f"GSD = {base}^({exponent_term_n} + {exponent_term_m} - {abs(exponent_constant)})")

solve_toric_code_gsd()