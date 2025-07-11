import math

def solve_arc_length_problem():
    """
    This function solves the given problem by following the plan outlined above.
    It calculates the length of the astroid's arc in the first quadrant and,
    based on that, determines the value of 'a'.
    """
    
    # Step 1-4: The integrand for arc length is |3*cos(t)*sin(t)|.
    # For the first quadrant (0 <= t <= pi/2), this is 3*cos(t)*sin(t).
    
    # Step 5: Calculate the length of the arc in the first quadrant.
    # We can integrate 3*sin(t)*cos(t) from 0 to pi/2.
    # The antiderivative of 3*sin(t)*cos(t) is (3/2)*sin(t)^2.
    
    t1 = 0
    t2 = math.pi / 2
    
    length_t1 = (3/2) * (math.sin(t1)**2)
    length_t2 = (3/2) * (math.sin(t2)**2)
    
    arc_length = length_t2 - length_t1
    
    print(f"The length of the arc from t={t1} to t={t2:.4f} is calculated as:")
    print(f"L = (3/2)*sin({t2:.4f})^2 - (3/2)*sin({t1})^2")
    print(f"L = {length_t2:.2f} - {length_t1:.2f} = {arc_length:.2f}")
    
    # Given arc length in the problem
    given_length = 3/2
    
    print(f"\nThis calculated length {arc_length:.2f} matches the given length {given_length:.2f}.")
    print("Therefore, we can conclude the arc is the portion in the first quadrant, corresponding to t from 0 to pi/2.")
    
    # Step 6: Determine the range of x-values for this arc.
    # x(t) = cos(t)^3
    x_at_t1 = math.cos(t1)**3
    x_at_t2 = math.cos(t2)**3
    
    # The range is from the minimum to the maximum of these values.
    x_min = min(x_at_t1, x_at_t2)
    x_max = max(x_at_t1, x_at_t2)
    
    print("\nFor this arc, the x-values range from x_min to x_max.")
    print(f"x_max corresponds to t={t1}: cos({t1})^3 = {int(x_at_t1)}")
    print(f"x_min corresponds to t={t2:.4f}: cos({t2:.4f})^3 = {int(x_at_t2)}")
    print(f"So the range of x is [{int(x_min)}, {int(x_max)}].")
    
    # Step 7: The problem states the range is [0, a].
    # By comparing [0, 1] with [0, a], we get a = 1.
    a = x_max
    
    print("\nThe problem defines this range as [0, a].")
    print("Therefore, the final equation is:")
    print(f"a = {int(a)}")
    print("\nThe value of 'a' is 1.")

solve_arc_length_problem()