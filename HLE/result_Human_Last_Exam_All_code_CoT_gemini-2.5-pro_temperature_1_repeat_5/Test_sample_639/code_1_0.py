import math

def calculate_packing_distances():
    """
    Calculates the possible normalized distances between hard spheres on a 2D plane for r <= 3.

    These distances correspond to peaks in the radial distribution function g(r)
    and are derived from the geometry of square and hexagonal packings.
    """
    # The set of unique squared normalized distances (r'^2) for r' <= 3 is found by
    # combining the possible distances from both square and hexagonal lattices.
    # r'^2 = n^2 + m^2 (square lattice)
    # r'^2 = n^2 + m^2 + nm (hexagonal lattice)
    #
    # Examples:
    # r'^2 = 1: (1^2 + 0^2) contact distance
    # r'^2 = 2: (1^2 + 1^2) square diagonal
    # r'^2 = 3: (1^2 + 1^2 + 1*1) hexagonal long diagonal
    # r'^2 = 4: (2^2 + 0^2) two diameters
    # r'^2 = 5: (2^2 + 1^2) knight's move on square lattice
    # r'^2 = 7: (2^2 + 1^2 + 2*1) from hexagonal lattice
    # r'^2 = 8: (2^2 + 2^2) diagonal of 2x2 square
    # r'^2 = 9: (3^2 + 0^2) three diameters
    
    squared_distances = [1, 2, 3, 4, 5, 7, 8, 9]
    
    distances = []
    print("The set of normalized distances r' = r/d for r' <= 3 are:")
    
    for r_squared in sorted(squared_distances):
        r = math.sqrt(r_squared)
        distances.append(r)
        # The prompt requested to "output each number in the final equation"
        # The following line shows the calculation for each distance.
        print(f"sqrt({r_squared}) = {r:.2f}")

if __name__ == "__main__":
    calculate_packing_distances()