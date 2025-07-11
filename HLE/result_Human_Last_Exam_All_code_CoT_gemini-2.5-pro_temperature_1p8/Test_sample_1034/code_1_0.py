import math

def solve_triangle_problem():
    """
    Calculates the number of integer values for a side of a triangle in the SSA ambiguous case.

    Let S = SNYWR. In triangle ABC, BC = S and m∠B = 30°. 
    This function computes the number of integer values of AC for which there are 
    two possible values for side length AB.
    """
    
    # Step 1: Set the value for S from SNYWR.
    S = 49197
    side_BC = S
    angle_B_deg = 30

    print(f"Given triangle ABC with BC = S = {S} and angle B = {angle_B_deg} degrees.")
    print("We need to find the number of integer values for side AC that result in two possible lengths for side AB.")
    print("-" * 30)
    
    # Step 2: Formulate the condition for the ambiguous case (SSA).
    # For two possible triangles, the length of the side opposite the given angle (AC)
    # must be greater than the altitude (h) from vertex C to side AB, and less than
    # the adjacent side (BC).
    # Condition: h < AC < BC
    # Where h = BC * sin(B)
    
    print("This is the ambiguous SSA case. For two solutions to exist, the length of side AC must be")
    print("strictly between the altitude from vertex C and the length of side BC.")
    
    # Step 3: Calculate the bounds for the length of AC.
    # sin(30 degrees) = 0.5
    h = side_BC * 0.5
    lower_bound = h
    upper_bound = side_BC

    print(f"\nThe required inequality for AC is: BC * sin(B) < AC < BC")
    print(f"Calculating the bounds: {side_BC} * sin({angle_B_deg}°) < AC < {side_BC}")
    print(f"{side_BC * 0.5} < AC < {side_BC}")
    print("-" * 30)

    # Step 4: Count the number of integers within these bounds.
    # The integers must be in the range (24598.5, 49197).
    
    # The first integer greater than the lower bound
    first_integer = math.floor(lower_bound) + 1
    
    # The last integer less than the upper bound
    last_integer = upper_bound - 1
    
    # Calculate the total count
    count = last_integer - first_integer + 1
    
    print("We now count the number of integers AC in this range.")
    print(f"The smallest possible integer for AC is {first_integer}.")
    print(f"The largest possible integer for AC is {last_integer}.")
    
    # Step 5: Display the final calculation as an equation.
    print("\nThe number of integer values is 'Last Integer - First Integer + 1'.")
    print(f"Final calculation: {last_integer} - {first_integer} + 1 = {count}")

solve_triangle_problem()