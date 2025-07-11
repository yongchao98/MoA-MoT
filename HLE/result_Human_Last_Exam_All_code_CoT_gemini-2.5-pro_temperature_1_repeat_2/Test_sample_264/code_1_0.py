import fractions

def solve_intersection_volume():
    """
    Calculates the volume of the intersection of two shifted unit cubes.
    """
    # Side length of the original cubes
    side_length = 1.0

    # The shift along each axis
    shift = 0.5

    # Define the boundaries of the first cube (C1)
    c1_min = 0.0
    c1_max = side_length

    # Define the boundaries of the second, shifted cube (C2)
    c2_min = shift
    c2_max = side_length + shift

    # The intersection is a region where a point is in both cubes.
    # We find the boundaries of this intersection for one dimension (e.g., x-axis).
    # The logic is the same for the y and z axes.
    intersect_min = max(c1_min, c2_min)
    intersect_max = min(c1_max, c2_max)

    # Calculate the side length of the resulting intersection cube
    intersect_side_length = intersect_max - intersect_min

    # Calculate the volume of the intersection
    volume = intersect_side_length ** 3

    # Print the explanation and the final result
    print("The first cube occupies the space from 0 to 1 on each axis.")
    print(f"The second cube is shifted by {shift} and occupies the space from {shift} to {side_length + shift} on each axis.")
    print("\nThe intersection of these two cubes is a smaller cube.")
    print(f"The side length of this intersection cube is the overlap on any axis, which is {intersect_max} - {intersect_min} = {intersect_side_length}.")
    print("\nThe volume is calculated by cubing this side length.")
    
    # Using fractions for a clear representation of the equation
    s_frac = fractions.Fraction(intersect_side_length)
    v_frac = fractions.Fraction(volume)

    print(f"\nFinal Equation: ({s_frac.numerator}/{s_frac.denominator}) * ({s_frac.numerator}/{s_frac.denominator}) * ({s_frac.numerator}/{s_frac.denominator}) = {v_frac.numerator}/{v_frac.denominator}")
    print(f"\nThe volume of the intersection is {float(volume)}.")

solve_intersection_volume()