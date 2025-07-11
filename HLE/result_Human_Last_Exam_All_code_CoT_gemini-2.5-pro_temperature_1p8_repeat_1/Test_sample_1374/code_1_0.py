import math

def solve_playdough_gravity():
    """
    Calculates the furthest distance from point A for a playdough shape
    that maximizes gravity at A.
    The derivation leads to the formula: r_max = (15 / (4 * pi))^(1/3).
    """
    numerator = 15
    denominator_const = 4
    pi = math.pi

    print("To find the furthest point on the surface of the playdough, we use the derived formula:")
    print(f"r_max = (numerator / (denominator_const * π))^(1/3)")
    print("\n--- Calculation Steps ---")
    
    # Print the values used in the equation
    print(f"1. Numerator: {numerator}")
    print(f"2. Denominator constant: {denominator_const}")
    print(f"3. Value of π: {pi}")

    # Calculate the term inside the parenthesis
    term_inside_power = numerator / (denominator_const * pi)
    print(f"4. The value of ({numerator} / ({denominator_const} * π)) is: {term_inside_power}")

    # Calculate the final distance
    final_distance = term_inside_power**(1/3)
    print(f"5. The final distance is the cube root of the above value.")
    
    print("\n--- Final Answer ---")
    print(f"The furthest point on the surface of the playdough is {final_distance:.5f} meters from point A.")

solve_playdough_gravity()