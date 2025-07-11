import sympy
from sympy import cos, sqrt, acos

def solve_billiards_limit():
    """
    This function calculates the limit of M(theta) as theta -> 0 for the given
    billiards problem using the sympy library.
    """
    # Define the symbolic variable theta
    theta = sympy.Symbol('theta')

    # As derived in the plan, the angle M(theta) is given by arccos of the
    # following expression, which represents the minimum possible value of cos(alpha).
    # This minimum is achieved when the distance between the start point x and
    # the impact point y is maximized (e.g., x=1, y=5*e^(i*theta)).
    cos_M_theta = (4 * cos(theta/2)) / sqrt(26 - 10 * cos(theta))

    # We want to find the limit of M(theta) = acos(cos_M_theta) as theta -> 0.
    # Since acos is continuous, we can compute the limit of the argument first.
    
    # Define the numerator and the term inside the square root in the denominator
    numerator = 4 * cos(theta/2)
    denominator_sq_term = 26 - 10 * cos(theta)
    
    # Calculate the limit of the numerator as theta -> 0
    limit_numerator = sympy.limit(numerator, theta, 0)
    
    # Calculate the limit of the term inside the square root as theta -> 0
    limit_denominator_sq_term = sympy.limit(denominator_sq_term, theta, 0)
    
    # Calculate the limit of the full expression for cos(M(theta))
    limit_cos_M = sympy.limit(cos_M_theta, theta, 0)

    # The final result is the arccosine of this limit
    final_limit = acos(limit_cos_M)

    # Print the step-by-step evaluation of the final equation
    print("The limit of M(theta) as theta approaches 0 is calculated as follows:")
    print(f"lim M(theta) = arccos( lim ( (4*cos(theta/2)) / sqrt(26 - 10*cos(theta)) ) )")
    print(f"             = arccos( (lim (4*cos(theta/2))) / sqrt(lim (26 - 10*cos(theta))) )")
    print(f"             = arccos( {limit_numerator} / sqrt({limit_denominator_sq_term}) )")
    print(f"             = arccos( {limit_numerator} / {sqrt(limit_denominator_sq_term)} )")
    print(f"             = arccos({limit_cos_M})")
    print(f"             = {final_limit}")

solve_billiards_limit()