import math

def find_planar_distances():
    """
    Calculates the possible normalized distances (r=R/d) between hard spheres
    in a 2D triangular lattice for r <= 3.

    In a 2D triangular lattice, the squared distance from the origin to any
    lattice point (n, m) is given by r^2 = n^2 + nm + m^2. We iterate
    through integer values of n and m to find all unique distances up to r=3.
    """
    r_max = 3
    # Use a set to store unique distances to avoid duplicates
    distances = set()

    # Iterate through a reasonable range of integer lattice coordinates (n, m)
    # Since r^2 = n^2 + nm + m^2, if n or m > 3, r is likely > 3.
    # A loop from 0 to 3 for both n and m is sufficient.
    for n in range(r_max + 1):
        for m in range(r_max + 1):
            # The formula for squared distance in a triangular lattice
            r_squared = float(n**2 + n * m + m**2)

            # Check if the distance is within the specified range (0 < r <= 3)
            # which is equivalent to (0 < r^2 <= 9)
            if 0 < r_squared <= r_max**2:
                distance = math.sqrt(r_squared)
                distances.add(distance)

    # Sort the unique distances in ascending order
    sorted_distances = sorted(list(distances))

    # Print the results formatted to two decimal places
    print("The set of possible normalized distances r for r <= 3 are:")
    # Create a list of formatted strings for the final answer
    result_list = []
    for r in sorted_distances:
        # We don't use the equation itself in the final printout,
        # but print the final numerical result as requested.
        print(f"{r:.2f}")
        result_list.append(f"{r:.2f}")

    # Final answer in the required format
    final_answer = ", ".join(result_list)
    print(f"\n<<<{final_answer}>>>")


if __name__ == "__main__":
    find_planar_distances()
