import math
import random

def solve_and_print():
    """
    This function calculates the area of the region R defined by the equation
    floor(|z|) = |floor(z)| within the 6x6 square in the complex plane.
    It does so by decomposing the total area into contributions from different
    unit squares and calculating each contribution using a Monte Carlo simulation.
    Finally, it prints the sum of these contributions and the total area.
    """
    
    # Properties of the unit squares that contribute to the area.
    # Format: (n, m, k, multiplicity) where k=sqrt(n^2+m^2) and multiplicity
    # accounts for symmetry (e.g., A(1,0) and A(0,1)).
    canonical_pairs_properties = [
        (0, 0, 0, 1),
        (1, 0, 1, 2),
        (2, 0, 2, 2),
        (3, 0, 3, 2),
        (4, 0, 4, 2),
        (5, 0, 5, 2),
        (3, 4, 5, 2)
    ]
    
    num_samples_per_square = 3000000
    
    total_area = 0.0
    area_parts = []

    for n, m, k, multiplicity in canonical_pairs_properties:
        count = 0
        for _ in range(num_samples_per_square):
            # Sample point (a,b) within the unit square [n, n+1) x [m, m+1)
            a = random.uniform(n, n + 1)
            b = random.uniform(m, m + 1)
            
            # Check the condition: floor(sqrt(a^2+b^2)) == k
            lhs = math.floor(math.sqrt(a**2 + b**2))
            if lhs == k:
                count += 1
        
        # Area of the valid region within one unit square [n,n+1)x[m,m+1).
        # The area of the square is 1, so the ratio is the area.
        sub_area = count / num_samples_per_square
        
        # Total contribution from this type of square (considering symmetry)
        total_contribution = multiplicity * sub_area
        area_parts.append(total_contribution)
        total_area += total_contribution
        
    # Print the final equation with each component rounded to two decimals
    equation_str = " + ".join([f"{part:.2f}" for part in area_parts])
    print(f"The total area is the sum of contributions from different unit squares:")
    print(f"Area = {equation_str} = {total_area:.2f}")

solve_and_print()