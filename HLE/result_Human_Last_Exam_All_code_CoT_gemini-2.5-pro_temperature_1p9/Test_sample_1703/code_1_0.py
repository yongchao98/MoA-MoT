import math

def get_area_formula():
    """
    This function explains the derivation of the area of the triangle T(t)
    and prints the final formula with its numerical components.
    """
    
    print("Step 1: The Area is Independent of Rotation")
    print("The rotation of the hexagon is a rigid body motion. It does not change the shape or the area of the inscribed triangle T(t). Therefore, we can calculate the area in the hexagon's own reference frame.\n")
    
    print("Step 2: The Triangle T(t) is Always Equilateral")
    print("At t=0, the vertices of T(t) are at the midpoints of alternating sides of the hexagon. The resulting triangle is equilateral with side length 15.")
    print("As the vertices move with the same speed along their respective sides, the symmetry is maintained, and the triangle T(t) remains equilateral. Its area only depends on its side length a(t).\n")
    
    print("Step 3: Finding the Area as a Function of Displacement")
    print("The area of an equilateral triangle is given by Area = (sqrt(3)/4) * a^2.")
    print("Through geometric analysis, the square of the triangle's side length, a(t)^2, can be expressed in terms of d(t), the displacement of the vertices from the midpoints of the hexagon's sides:")
    print("a(t)^2 = 225 + 3 * d(t)^2\n")
    
    print("Step 4: Modeling the Displacement d(t)")
    print("The vertices move with speed v=1 along the hexagon sides (length 10). Starting from the midpoint, a vertex reaches an endpoint (distance 5) in 5 seconds. The motion is a back-and-forth pattern.")
    print("The area depends on d(t)^2, which is periodic with a period of 10 seconds.")
    print("The function for the squared displacement is: d(t)^2 = (5 - |(t mod 10) - 5|)^2\n")
    
    print("Step 5: The Final Formula for the Area")
    print("Substituting the expression for d(t)^2 into the formula for a(t)^2, we get the final expression for the area as a function of time t:")
    print("Area(t) = (sqrt(3) / 4) * (225 + 3 * (5 - |(t mod 10) - 5|)^2)\n")

    print("--- Final Equation Breakdown ---")
    print("The final equation for the area is:")
    print("Area(t) = C1 * (C2 + C3 * (C4 - |(t mod C5) - C4|)^2)")
    print("The numerical values for the constants are:")
    
    constant_1_str = "sqrt(3) / 4"
    constant_2 = 225
    constant_3 = 3
    constant_4 = 5
    constant_5 = 10
    
    print(f"C1 = {constant_1_str}")
    print(f"C2 = {constant_2}")
    print(f"C3 = {constant_3}")
    print(f"C4 = {constant_4}")
    print(f"C5 = {constant_5}")

# Execute the function to print the solution
get_area_formula()

# The final formula in a compact form: Area(t) = (sqrt(3)/4) * (225 + 3*(5-abs((t%10)-5))**2)
# The value is Area(t) in square units.

# <<<Area(t) = (sqrt(3)/4) * (225 + 3 * (5 - abs((t % 10) - 5))**2)>>>