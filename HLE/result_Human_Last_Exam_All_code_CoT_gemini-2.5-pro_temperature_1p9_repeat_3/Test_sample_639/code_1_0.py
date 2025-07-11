import math

def calculate_planar_distances():
    """
    Calculates the set of possible normalized distances (r <= 3) between two hard
    spheres lying on the same plane in a jammed system.

    This considers arrangements corresponding to both hexagonal and square lattices
    to account for different possible local packings.
    """
    
    # Use a set to store unique distances
    distances = set()
    
    # Iterate through a grid of integer coefficients (n, m).
    # A range of -3 to 3 is sufficient since if n>3 or m>3, r will be > 3.
    limit = 3
    for n in range(-limit, limit + 1):
        for m in range(-limit, limit + 1):
            # Skip the origin point (distance to self)
            if n == 0 and m == 0:
                continue

            # 1. Calculate for hexagonal lattice
            # r^2 = n^2 + nm + m^2
            r_squared_hex = n**2 + n * m + m**2
            r_hex = math.sqrt(r_squared_hex)
            if r_hex <= 3:
                distances.add(r_hex)

            # 2. Calculate for square lattice
            # r^2 = n^2 + m^2
            r_squared_sq = n**2 + m**2
            r_sq = math.sqrt(r_squared_sq)
            if r_sq <= 3:
                distances.add(r_sq)

    # Sort the distances and format them
    sorted_distances = sorted(list(distances))
    
    # Format the numbers for the final equation string
    formatted_distances = [f"{d:.2f}" for d in sorted_distances]

    # Print the final equation as requested
    print("The set of possible normalized distances r for r <= 3 is:")
    # The problem asks to output each number in the final equation.
    # The format "r = {d1, d2, ...}" seems appropriate.
    print(f"r = {{{', '.join(formatted_distances)}}}")

if __name__ == "__main__":
    calculate_planar_distances()