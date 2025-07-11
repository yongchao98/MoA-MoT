import math

def calculate_limit():
    """
    This script calculates the limit of M(theta) as theta approaches 0.
    
    The analysis shows that the supremum of the angle alpha, M(theta),
    is achieved when cos(alpha) is minimized. For theta -> 0, this occurs
    when the trajectory starts near one endpoint of the allowed arc |x|=1
    and ends at the farthest vertex of the side A.

    Let's take the starting point x to be (1, 0) and the endpoint y
    to be the vertex Q = (5*cos(theta), 5*sin(theta)).
    
    The expression for cos(M(theta)) in this limiting case is:
    cos(M(theta)) = (4 * cos(theta/2)) / sqrt(26 - 10*cos(theta))

    We need to compute the limit of M(theta) as theta -> 0, which is:
    lim_{theta->0} M(theta) = arccos( lim_{theta->0} cos(M(theta)) )
    
    The script below performs this final calculation step-by-step.
    """

    print("We want to find the limit of M(theta) as theta goes to 0.")
    print("M(0) = lim_{theta->0} arccos((4*cos(theta/2)) / sqrt(26 - 10*cos(theta)))")
    print("This can be solved by finding the limit of the argument first, then taking the arccos.")

    print("\nStep 1: Evaluate the limit of the argument of arccos.")
    print("lim_{theta->0} (4*cos(theta/2)) / sqrt(26 - 10*cos(theta))")
    
    # We substitute theta = 0 into the expression.
    cos_0 = math.cos(0)
    print(f"= (4 * cos(0)) / sqrt(26 - 10*cos(0))")
    print(f"= (4 * {cos_0}) / sqrt(26 - 10*{cos_0})")
    
    numerator = 4 * cos_0
    denominator_radicand = 26 - 10 * cos_0
    print(f"= {numerator} / sqrt({denominator_radicand})")
    
    denominator = math.sqrt(denominator_radicand)
    print(f"= {numerator} / {denominator}")
    
    limit_value = numerator / denominator
    print(f"= {limit_value}")

    print("\nStep 2: Evaluate arccos of the limit.")
    
    # Take the arccosine of the result
    final_result_rad = math.acos(limit_value)
    print(f"M(0) = arccos({limit_value})")
    print(f"M(0) = {final_result_rad} radians")
    
# Execute the calculation
calculate_limit()

# Final answer in the required format
print("<<<0>>>")