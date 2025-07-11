import math

def solve_triangle_area():
    """
    This function derives and prints the formula for the area of triangle T(t).
    """
    # Given parameters
    # Radius of the circumscribed circle C
    R = 10.0
    # Speed of the vertices of T(t)
    v = 1.0
    
    # 1. The rotation of the hexagon (omega) is irrelevant for the area calculation,
    # as area is invariant under rotation.
    
    # 2. Due to symmetry, the triangle T(t) is always equilateral.
    # The area is Area = (sqrt(3)/4) * L^2, where L is the side length.
    
    # 3. We found the side length squared L(t)^2 is a function of d(t),
    # the displacement from the midpoint of the hexagon side:
    # L(t)^2 = (9*R^2/4) + 3*d(t)^2
    
    # Calculate the constant term in the L(t)^2 expression
    # This corresponds to the area at t=0 when d(t)=0.
    const_term = (9 * R**2) / 4
    
    # The coefficient of the d(t)^2 term.
    d_squared_coeff = 3.0
    
    # 4. Model the displacement d(t). The motion is periodic.
    # A vertex moves from the midpoint to an end and back. The path length
    # for one cycle of displacement (e.g., mid -> end -> mid) is R.
    # The period of the displacement magnitude is T = R / v.
    # The amplitude is A = R / 2.
    # The function for d(t) is a triangular wave: d(t) = A - abs((t % T) - A)
    d_amplitude = R / 2
    d_period = R / v
    d_offset = R / (2 * v) # which is equal to d_amplitude in this case

    # 5. Assemble the final equation for the area.
    # We will print out each component of the final equation.
    
    print("The area of the triangle T(t) is given by the function:")
    print(f"Area(t) = (sqrt(3)/4) * ({const_term} + {d_squared_coeff} * d(t)^2)")
    print("\nwhere d(t) is the periodic displacement from the midpoint, given by:")
    print(f"d(t) = {d_amplitude} - abs((t % {d_period}) - {d_offset})")
    
    print("\nTherefore, the final combined equation for the area is:")
    # We use integers for a cleaner final look where possible.
    final_equation = f"Area(t) = (sqrt(3)/4) * ({int(const_term)} + {int(d_squared_coeff)} * ({int(d_amplitude)} - abs((t % {int(d_period)}) - {int(d_offset)}))^2)"
    print(final_equation)

solve_triangle_area()