import math

def calculate_area():
    """
    Calculates the area of the region R based on the analytical solution derived from its definition.
    The region R is where floor(|z|) = |floor(z)|. This equality holds only if |floor(z)| is an integer, k.
    This restricts the problem to unit squares (nx, ny) where nx^2+ny^2 is a perfect square.
    The condition for points z in such a square becomes k <= |z| < k+1.
    The total area is the sum of areas of intersections of these annuli with their corresponding unit squares.
    """
    pi = math.pi
    sqrt = math.sqrt
    asin = math.asin

    # List of terms for the total area. We use symmetry for pairs like (1,0) and (0,1).
    area_components = []

    # Case (0,0): k=0. Intersection of [0,1]x[0,1] with annulus 0 <= |z| < 1.
    # This is a quarter circle of radius 1. Area = pi * 1^2 / 4
    area_0_0 = pi / 4
    area_components.append(area_0_0)

    # Case (1,0) & (0,1): k=1. Annulus 1 <= |z| < 2.
    # Area per square: Integral from 1 to 2 of sqrt(4 - x^2) dx.
    area_1_0 = (2 * pi / 3 - sqrt(3) / 2)
    area_components.append(2 * area_1_0)
    
    # Case (2,0) & (0,2): k=2. Annulus 2 <= |z| < 3.
    # Area per square: Integral from 2 to 3 of sqrt(9 - x^2) dx.
    area_2_0 = (9 * pi / 4 - sqrt(5) - 4.5 * asin(2 / 3))
    area_components.append(2 * area_2_0)

    # Case (3,0) & (0,3): k=3. Annulus 3 <= |z| < 4.
    # Area per square: Integral from 3 to 4 of sqrt(16 - x^2) dx.
    area_3_0 = (4 * pi - 1.5 * sqrt(7) - 8 * asin(3 / 4))
    area_components.append(2 * area_3_0)

    # Case (4,0) & (0,4): k=4. Annulus 4 <= |z| < 5.
    # Area per square: Integral from 4 to 5 of sqrt(25 - x^2) dx.
    area_4_0 = (25 * pi / 4 - 6 - 12.5 * asin(4 / 5))
    area_components.append(2 * area_4_0)

    # Case (5,0) & (0,5): k=5. Annulus 5 <= |z| < 6.
    # Area per square: Integral from 5 to 6 of sqrt(36 - x^2) dx.
    area_5_0 = (9 * pi - 2.5 * sqrt(11) - 18 * asin(5 / 6))
    area_components.append(2 * area_5_0)

    # Case (3,4) & (4,3): k=5. Annulus 5 <= |z| < 6.
    # Area per square: Integral from 3 to 4 of (sqrt(36-x^2)-4) dx.
    area_3_4 = (4 * sqrt(5) - 4.5 * sqrt(3) - 4 + 18 * asin(2 / 3) - 3 * pi)
    area_components.append(2 * area_3_4)
    
    total_area = sum(area_components)
    
    # Output the final equation with each numerical component, as requested.
    equation_str = " + ".join([f"{v:.4f}" for v in area_components])
    print(f"Total Area = {equation_str}")
    print(f"Total Area = {total_area:.4f}")
    print(f"\nThe area of R, expressed to two decimals, is {total_area:.2f}")

calculate_area()