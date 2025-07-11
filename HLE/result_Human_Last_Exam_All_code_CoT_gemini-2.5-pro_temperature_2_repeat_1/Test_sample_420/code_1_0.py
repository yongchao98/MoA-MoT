import sympy

def solve_cube_area_problem():
    """
    Solves the cube surface area problem symbolically.

    Calculates the ratio of the area of region D to the total surface area S of a cube.
    D is the set of points on the surface S at a surface distance of at most sqrt(2)*s from a vertex P.
    """
    
    # Let 's' be the side length of the cube.
    # We can set s=1 for the ratio calculation as it will cancel out.
    s = 1

    # Define the mathematical constants we will need
    pi = sympy.pi
    sqrt3 = sympy.sqrt(3)

    # 1. Total surface area of the cube
    area_S = 6 * s**2
    print(f"Step 1: The total surface area of the cube, S, is 6 * s^2.")
    if s == 1:
        print(f"         For s=1, S = {area_S}\n")


    # 2. Area on the 3 faces adjacent to vertex P
    # The largest surface distance from P on an adjacent face is to the diagonally opposite
    # vertex of that face, which is sqrt(s^2 + s^2) = sqrt(2)*s.
    # This means the entire area of all three adjacent faces is included in D.
    area_D_adjacent = 3 * s**2
    print(f"Step 2: The distance to any point on the three faces adjacent to P is at most sqrt(2)*s.")
    print(f"         Thus, their entire area is included in D.")
    print(f"         Area on adjacent faces = 3 * s^2 = {area_D_adjacent} (for s=1)\n")


    # 3. Area on one of the 3 faces opposite to vertex P
    # The calculation for one opposite face involves finding the area on a square [0,s]x[0,s]
    # where min((s+y)^2+x^2, (s+x)^2+y^2) <= 2*s^2.
    # This is found using integral calculus and results in the following exact expression:
    area_D_one_opposite_face = s**2 * (2*pi - 9 + 3*sqrt3) / 6
    print(f"Step 3: The area on any single face opposite to P requires a more complex calculation involving")
    print(f"         unfolding the cube. The resulting area on one opposite face is:")
    print(f"         s^2 * (2*pi - 9 + 3*sqrt(3)) / 6\n")

    
    # 4. Total area of region D
    # Area(D) is the sum of the areas from the adjacent faces and the three opposite faces.
    area_D = area_D_adjacent + 3 * area_D_one_opposite_face
    print(f"Step 4: The total area of D is the sum of the area from the 3 adjacent faces")
    print(f"         and the 3 opposite faces.")
    print(f"         Area(D) = 3*s^2 + 3 * [s^2 * (2*pi - 9 + 3*sqrt(3)) / 6]")
    # Simplify the expression for Area(D)
    area_D_simplified = sympy.simplify(area_D)
    print(f"         Simplified, Area(D) = {area_D_simplified} (for s=1)\n")
    

    # 5. Calculate the final ratio
    ratio = area_D / area_S
    ratio_simplified = sympy.simplify(ratio)
    
    print(f"Step 5: The final ratio is Area(D) / Area(S).")
    print(f"         Ratio = ({area_D_simplified}) / ({area_S})")
    print(f"         The simplified ratio is: {ratio_simplified}\n")
    
    # Final step is to print the required output as an equation string
    # numerator = 2*pi - 3 + 3*sqrt(3), denominator = 12
    num, den = ratio_simplified.as_numer_denom()
    
    print("Final Answer Equation:")
    # We format the string to show the expression clearly
    print(f"({2} * pi - {3} + {3} * sqrt(3)) / {12}")

solve_cube_area_problem()