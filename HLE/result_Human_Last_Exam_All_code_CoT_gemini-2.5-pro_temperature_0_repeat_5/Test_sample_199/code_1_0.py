import math

def solve_minimal_polynomial():
    """
    Calculates the minimal polynomial for the shortest geodesic distance
    an ant can walk starting and ending at a vertex of a dodecahedron.
    """
    # Step 1: Define the golden ratio, phi.
    # phi = (1 + sqrt(5)) / 2
    
    # Step 2: Use the known result for the squared length (d^2) of the shortest
    # vertex-loop geodesic on a dodecahedron of side length 1.
    # d^2 = 12 * phi + 8
    # We can express this in the form a + b*sqrt(5).
    # d^2 = 12 * (1 + sqrt(5))/2 + 8 = 6 * (1 + sqrt(5)) + 8 = 14 + 6*sqrt(5)
    
    # Step 3: Simplify d = sqrt(14 + 6*sqrt(5)).
    # We look for d = a + b*sqrt(5).
    # d^2 = (a + b*sqrt(5))^2 = a^2 + 5b^2 + 2ab*sqrt(5).
    # Comparing with 14 + 6*sqrt(5), we have:
    # 2ab = 6  => ab = 3
    # a^2 + 5b^2 = 14
    # By inspection (or solving), a=3 and b=1 are integer solutions.
    # So, d = 3 + sqrt(5).
    a = 3
    b = 1
    
    # Step 4: Find the minimal polynomial for d = a + b*sqrt(5).
    # Let x = d.
    # x - a = b*sqrt(5)
    # (x - a)^2 = (b^2)*5
    # x^2 - 2ax + a^2 - 5b^2 = 0
    
    coeff_x2 = 1
    coeff_x1 = -2 * a
    coeff_const = a**2 - 5 * b**2
    
    # Step 5: Print the final equation with each coefficient.
    print("The minimal polynomial is:")
    print(f"{coeff_x2}x^2 + ({coeff_x1})x + {coeff_const} = 0")

solve_minimal_polynomial()