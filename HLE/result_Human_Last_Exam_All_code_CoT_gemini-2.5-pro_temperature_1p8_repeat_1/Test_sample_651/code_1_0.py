import sympy
from sympy import Symbol, pi, limit

def solve_limit():
    """
    This function outlines the symbolic solution to the problem.

    The derivation involves the following key findings for small theta > 0:
    1. The range of possible trajectory angles `phi` is approximately [-theta/4, 5*theta/4].
       This is found by analyzing the vectors connecting the endpoints of the starting arc
       (from 1 to e^(i*theta)) and the target segment A (from 5 to 5*e^(i*theta)).
    2. The direction of the inner normal vector `psi` is approximately pi - theta/2.
    3. The angle `alpha` between the trajectory and the inner normal is |phi - psi|.
       Physical constraints require the angle to be obtuse.
    4. The set of all possible alpha values is the interval [pi - (7/4)*theta, pi - (1/4)*theta].
    5. The supremum M(theta) is the largest possible value for alpha.
    """
    
    # Define theta as a symbolic variable
    theta = Symbol('theta', real=True, positive=True)

    # For small theta, M(theta) is the maximum possible angle alpha.
    # From our derivation, this corresponds to the trajectory from x = exp(i*theta) to y = 5.
    # The angle alpha for this trajectory is |arg(5 - exp(i*theta)) - arg(inner_normal)|.
    # Approximating this for small theta gives:
    # M(theta) approx |(-theta/4) - (pi - theta/2)| = |-pi + theta/4| = pi - theta/4.
    
    num = 1
    den = 4
    M_theta = pi - (num / den) * theta

    # Calculate the limit of M(theta) as theta approaches 0.
    final_limit = limit(M_theta, theta, 0)
    
    # Output the explanation and the result.
    print("The supremum of the angle alpha, M(theta), for small positive theta is given by the expression:")
    # The prompt requests printing numbers in the equation.
    print(f"M(theta) = pi - ({num}/{den}) * theta")
    print("\nThe limit of M(theta) as theta approaches 0 is:")
    print(f"lim_{{theta->0}} (pi - {num}*theta/{den}) = {final_limit}")

solve_limit()
