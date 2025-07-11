import math

def solve_triangle_problem():
    """
    Calculates the number of integer values for side AC for which there are 
    two possible triangles ABC, given BC and angle B.
    """
    
    # Step 1: Decode the value of S from SNYWR
    # S=3, N=6, Y=6, W=2, R=4
    s_value = 36624
    
    # Step 2: Assign the given parameters
    bc_length = s_value
    angle_b_degrees = 30
    
    # Step 3: Formulate the inequality for the ambiguous case (SSA)
    # The condition for two possible triangles is: BC * sin(B) < AC < BC
    
    # Step 4: Calculate the lower and upper bounds for AC's length
    # sin(30 degrees) is exactly 0.5
    lower_bound = bc_length * math.sin(math.radians(angle_b_degrees))
    upper_bound = bc_length
    
    # We need to count the number of integers AC such that:
    # lower_bound < AC < upper_bound
    # Note: Since lower_bound is an integer, the smallest integer value for AC
    # is lower_bound + 1. The largest is upper_bound - 1.
    
    # Step 5: Count the integers in the calculated range.
    # The number of integers between x and y (exclusive) is y - x - 1.
    # We use int() to ensure we are working with integer boundaries.
    num_integer_values = int(upper_bound) - int(lower_bound) - 1
    
    # Output the explanation and the final equation
    print(f"The value of S (SNYWR) is {s_value}.")
    print(f"The length of side BC is {bc_length}.")
    print(f"The measure of angle B is {angle_b_degrees} degrees.")
    print("\nFor there to be two possible triangles, the length of side AC must be in the interval (BC * sin(B), BC).")
    print(f"The lower bound for AC is BC * sin({angle_b_degrees}) = {bc_length} * 0.5 = {lower_bound}.")
    print(f"The upper bound for AC is {upper_bound}.")
    print(f"So, the length of AC must satisfy: {int(lower_bound)} < AC < {int(upper_bound)}.")
    
    print("\nThe number of integers in this range is found by the calculation:")
    print(f"Number of values = Upper Bound - Lower Bound - 1")
    print(f"Number of values = {int(upper_bound)} - {int(lower_bound)} - 1 = {num_integer_values}")

solve_triangle_problem()
