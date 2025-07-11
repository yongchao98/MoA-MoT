import math
# The fpylll library is required for this solution.
# Please install it using: pip install fpylll
# mpmath is also required: pip install mpmath
from fpylll import IntegerMatrix, LLL
from mpmath import mp

def solve_arctan_relation():
    """
    Finds the integer coefficients for the given Machin-like formula
    using the LLL algorithm.
    """
    # Set precision for floating point operations to 150 digits.
    mp.dps = 150

    # Define the input values from the problem statement.
    xs = [122, 239, 682, 1252, 2855, 12943]

    # We are searching for an integer relation for the vector:
    # [arctan(1/122), ..., arctan(1/12943), -pi/4]
    reals = [mp.atan(mp.mpf(1)/x) for x in xs]
    reals.append(-mp.pi / 4)

    dim = len(reals)

    # A large integer constant for scaling the real values into integers.
    C = mp.power(10, mp.dps - 5)

    # Create the basis matrix for the LLL algorithm to find the integer relation.
    # The lattice is constructed such that a short vector in its basis
    # corresponds to the integer coefficients of the linear relation.
    B = IntegerMatrix(dim, dim + 1)
    for i in range(dim):
        B[i, i] = 1

    for i in range(dim):
        B[i, dim] = int(C * reals[i])

    # Apply the LLL algorithm to find a reduced basis for the lattice.
    B_reduced = LLL.reduction(B)

    # The first row of the reduced basis is the shortest vector found by LLL.
    # The first `dim` components of this vector are the integer coefficients.
    relation_coeffs = list(B_reduced[0])[:dim]

    # The coefficients correspond to [c1, c2, c3, c4, c5, c6, n_temp].
    # The equation is: c1*atan(...) + ... + n_temp*(-pi/4) = 0
    # which is equivalent to: c1*atan(...) + ... = n_temp*(pi/4)
    # So, the desired n is n_temp.
    n = relation_coeffs[-1]
    constants = relation_coeffs[:-1]

    # The problem requires the smallest positive n.
    # If LLL returns a relation with a negative n, we can flip the signs
    # of all coefficients to get an equivalent relation with a positive n.
    if n < 0:
        n = -n
        constants = [-c for c in constants]
    elif n == 0:
        # This is unlikely for a well-posed problem, but as a safeguard:
        print("Warning: Found a solution with n=0. This might not be the desired one.")


    # Print the full equation as requested.
    c1, c2, c3, c4, c5, c6 = constants
    
    print("The discovered relation is:")
    
    equation_str = f"{n} * pi/4 ="
    
    # Add terms for each constant and its corresponding arctan
    # The logic here handles positive and negative coefficients to build the string.
    is_first_term = True
    for i, c in enumerate(constants):
        if c == 0:
            continue
        
        # Add sign between terms
        if not is_first_term:
            if c > 0:
                equation_str += " + "
            else:
                equation_str += " - "
        else:
             if c < 0:
                equation_str += " -"

        is_first_term = False
        
        # Add the coefficient and arctan term
        equation_str += f" {abs(c)}*arctan(1/{xs[i]})"

    print(equation_str)
    print()

    # Print the final answer in the requested format
    final_answer = [n] + constants
    print("The solution for (n, c1, c2, c3, c4, c5, c6) is:")
    print(','.join(map(str, final_answer)))

if __name__ == '__main__':
    solve_arctan_relation()
