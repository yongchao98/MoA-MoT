import math

def solve():
    """
    This function explains and prints the condition for the solution to blow up.
    The problem asks for the values of y(0) for which the solution blows up,
    given x(0) > 1. The analysis leads to a condition that depends on the
    initial value x(0).

    The condition is determined by a separatrix in the phase plane, which is
    the stable manifold of the saddle point at (1,0). The equation for this
    separatrix is y^2 = 2*x + 1 - 3*x^(2/3).

    For x > 1, the solution blows up if the initial point (x(0), y(0)) lies
    below the upper branch of this separatrix.
    """

    print("For a given initial value x(0) > 1, the solution to the system of differential equations blows up if y(0) satisfies the following inequality:")

    # The final equation is y(0) < (x(0)**(1/3) - 1) * sqrt(2*x(0)**(1/3) + 1)
    # The numbers in this equation are:
    power_val_numerator = 1
    power_val_denominator = 3
    term1_const = -1
    term2_coeff = 2
    term2_const = 1

    # To meet the requirement of outputting each number in the final equation,
    # we construct the expression from these numbers.
    # We use x_0 as a placeholder for x(0) in the printout for clarity.
    
    power_str = f"{power_val_numerator}/{power_val_denominator}"

    # Construct and print the final inequality.
    print("\nIn symbolic form:")
    print(f"y(0) < (x(0)^({power_str}) + ({term1_const})) * sqrt({term2_coeff}*x(0)^({power_str}) + {term2_const})")

solve()
