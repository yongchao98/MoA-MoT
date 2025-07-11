import math

def solve_triangle_area():
    """
    Calculates and prints the formula for the area of the triangle T(t).
    """

    # Step 1: Define the properties of the hexagon.
    # A regular hexagon is inscribed in a circle of radius R.
    # Its side length 's' is equal to the radius.
    R = 10
    s = R

    # Step 2: Determine the properties of the initial triangle T(0).
    # Its vertices are midpoints of alternating sides (e.g., AB, CD, EF).
    # We can calculate the side length L(0) of this triangle. By placing the
    # hexagon's center at the origin, we can find the coordinates of the midpoints.
    # Midpoint of side 1 (M1): (3*R/4, R*math.sqrt(3)/4)
    # Midpoint of side 2 (M2): (-3*R/4, R*math.sqrt(3)/4)
    # The distance between M1 and M2 is the side length L(0).
    # L0 = (3*R/4) - (-3*R/4) = 6*R/4 = 1.5 * R
    L0 = 1.5 * R
    
    print("Step-by-step derivation of the area A(t):")
    print(f"1. The hexagon has side length s = R = {R}.")
    print(f"2. The initial triangle T(0) is equilateral with side length L(0) = 1.5 * R = {L0}.")
    
    # Step 3: Analyze the motion of the vertices.
    # The vertices move from the midpoints with speed v=1.
    # The distance moved from the midpoint at time t is d = v*t = t.
    v = 1
    
    # The squared side length of the triangle T(t) is L(t)^2 = L(0)^2 + 3*(v*t)^2
    L0_squared = L0**2
    t_squared_coefficient = 3 * (v**2)

    print(f"3. The vertices move a distance d=v*t={v}*t=t from the side midpoints.")
    print("4. The rotation of the hexagon does not change the triangle's area, as area is invariant under rotation.")
    print("5. The squared side length of the triangle at time t is given by L(t)^2.")
    print(f"   L(t)^2 = L(0)^2 + 3*d^2 = {L0}^2 + 3*t^2")
    print(f"   L(t)^2 = {L0_squared} + {t_squared_coefficient}*t^2")

    # Step 4: Calculate the area A(t).
    # The area of an equilateral triangle with side L is (sqrt(3)/4) * L^2.
    denominator = 4

    print("\n6. The area of an equilateral triangle is Area = (sqrt(3)/4) * L^2.")
    print("   Substituting L(t)^2 gives the final formula for A(t):")
    
    final_formula = f"A(t) = (sqrt(3) / {denominator}) * ({L0_squared} + {t_squared_coefficient}*t^2)"
    expanded_formula = f"A(t) = ({L0_squared}*sqrt(3) / {denominator}) + ({t_squared_coefficient}*sqrt(3) / {denominator})*t^2"
    
    print("\n" + "="*50)
    print("Final Formula:")
    print(final_formula)
    print("\nAlternatively, in expanded form:")
    print(expanded_formula)
    print("="*50)

solve_triangle_area()