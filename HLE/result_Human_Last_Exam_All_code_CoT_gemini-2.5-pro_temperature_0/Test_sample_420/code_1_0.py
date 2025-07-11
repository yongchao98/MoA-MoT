import math

def solve_cube_area_problem():
    """
    Calculates the ratio of the area of region D to the total surface area S of a cube.

    The solution is derived step-by-step:
    1. The area on the 3 faces adjacent to vertex P is calculated.
    2. The area on the 3 faces opposite to vertex P is calculated by considering
       the union of two circular segments resulting from unfolding the cube.
    3. The total area of D is summed up.
    4. The ratio Area(D) / Area(S) is computed and simplified to its exact form.
    """
    # Let s=1 for simplicity, as it will cancel out in the final ratio.
    s = 1
    
    # 1. Total surface area of the cube
    area_S_val = 6 * s**2
    
    # 2. Area on the 3 faces adjacent to P
    # As derived in the plan, the distance to any point on these faces is <= sqrt(2)*s.
    # So, the entire area of these 3 faces is included.
    area_adjacent_faces = 3 * s**2
    
    # 3. Area on one of the 3 opposite faces
    # This area is the union of two regions, R1 and R2.
    # Area(R1) = Area(R2) = s^2 * (pi/4 - 1/2)
    # Area(R1 intersect R2) can be calculated through integration.
    # A_int = (s**2 / 6) * (pi + 3 - 3*sqrt(3))
    # Area_opp = Area(R1) + Area(R2) - A_int
    # Area_opp = s**2 * (pi/2 - 1) - (s**2 / 6) * (pi + 3 - 3*sqrt(3))
    # Area_opp = (s**2 / 6) * (3*pi - 6 - pi - 3 + 3*sqrt(3))
    # Area_opp = (s**2 / 6) * (2*pi - 9 + 3*sqrt(3))
    
    # Let's define the components for the numerator of Area_opp
    # Numerator is (2*pi - 9 + 3*sqrt(3))
    # Denominator is 6
    
    # 4. Total area of D
    # Area(D) = area_adjacent_faces + 3 * Area_opp
    # Area(D) = 3*s^2 + 3 * (s**2 / 6) * (2*pi - 9 + 3*sqrt(3))
    # Area(D) = 3*s^2 + (s**2 / 2) * (2*pi - 9 + 3*sqrt(3))
    # Area(D) = s**2 * (3 + pi - 9/2 + (3*sqrt(3))/2)
    # Area(D) = s**2 * (pi - 3/2 + (3*sqrt(3))/2)
    # Area(D) = (s**2 / 2) * (2*pi - 3 + 3*sqrt(3))
    
    # 5. The final ratio Area(D) / Area(S)
    # Ratio = Area(D) / (6*s^2)
    # Ratio = [(s**2 / 2) * (2*pi - 3 + 3*sqrt(3))] / (6*s**2)
    # Ratio = (2*pi - 3 + 3*sqrt(3)) / 12
    
    # The components of the final fraction's numerator and denominator
    num_pi_coeff = 2
    num_const_term = -3
    num_sqrt3_coeff = 3
    denominator = 12
    
    print("The area of region D divided by the area of the surface S is given by the fraction:")
    print(f"  ( {num_pi_coeff}*pi + ({num_const_term}) + {num_sqrt3_coeff}*sqrt(3) )")
    print(f"--------------------------------------")
    print(f"                 {denominator}          ")

    # For verification, let's calculate the numerical value
    numerical_value = (num_pi_coeff * math.pi + num_const_term + num_sqrt3_coeff * math.sqrt(3)) / denominator
    print(f"\nNumerical approximation: {numerical_value}")

solve_cube_area_problem()
<<< (2*pi - 3 + 3*sqrt(3)) / 12 >>>