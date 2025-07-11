import math

def calculate_planar_distances():
    """
    Calculates the set of possible normalized distances (r <= 3) between
    hard spheres arranged on a 2D triangular lattice.
    """
    # Use a set to store unique distances to avoid duplicates.
    # The normalized distance r squared is given by r^2 = n^2 + n*m + m^2
    # for integers n, m on a triangular lattice. We are looking for r <= 3, so r^2 <= 9.
    distances = set()

    # The search range for n and m can be limited. If |n|>3 or |m|>3, r^2 will likely be > 9.
    # A safe range is -3 to 3 for both.
    n_range = range(-3, 4)
    m_range = range(-3, 4)

    for n in n_range:
        for m in m_range:
            # The case (n,m) = (0,0) is the distance from a particle to itself, which is 0.
            if n == 0 and m == 0:
                continue

            r_squared = n**2 + n * m + m**2

            if r_squared <= 9:
                r = math.sqrt(r_squared)
                distances.add(r)

    # Sort the distances and print them with the required formatting.
    sorted_distances = sorted(list(distances))
    
    print("The set of possible normalized distances r <= 3 are:")
    for r in sorted_distances:
        # The prompt requires "you still need to output each number in the final equation!".
        # We interpret this as just printing the final resulting numbers.
        print(f"{r:.2f}")

calculate_planar_distances()
