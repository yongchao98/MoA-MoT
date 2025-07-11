import math

def calculate_triangle_area_formula():
    """
    This script calculates the coefficients for the formula of the area of triangle T(t).
    """

    # --- Given Parameters ---
    # Radius of the circumscribed circle C
    R = 10
    # Speed of the triangle's vertices along the hexagon's sides
    v = 1

    # --- Step 1: Define Hexagon Properties ---
    # The side length 's' of a regular hexagon inscribed in a circle is equal to the radius R.
    s = R

    # --- Step 2: Derive the Squared Side Length of the Triangle T(t) ---
    # Based on geometric derivation, the squared side length L(t)^2 of the equilateral
    # triangle T(t) is given by the formula:
    # L(t)^2 = (9 * s^2 / 4) + 3 * (v*t)^2
    # We will now calculate the two main coefficients of this polynomial in t.
    
    # The constant term of L(t)^2
    constant_term = (9/4) * (s**2)

    # The coefficient of the t^2 term
    t_squared_coefficient = 3 * (v**2)

    # --- Step 3: Formulate the Area Expression ---
    # The area of an equilateral triangle is (sqrt(3)/4) * L^2.
    # We substitute our expression for L(t)^2.
    # Area(t) = (sqrt(3) / 4) * (constant_term + t_squared_coefficient * t^2)

    # The problem asks for the final equation with each number outputted.
    # We will print the final formula using the calculated coefficients.
    
    print("The area of the triangle T(t) as a function of time is given by the formula:")
    
    # Storing the numbers of the final equation
    num_sqrt_3 = "sqrt(3)"
    num_4 = 4
    num_const = int(constant_term)
    num_t_coeff = int(t_squared_coefficient)

    # Printing the final equation with the derived numbers
    print(f"Area(t) = ({num_sqrt_3} / {num_4}) * ({num_const} + {num_t_coeff} * t^2)")

calculate_triangle_area_formula()