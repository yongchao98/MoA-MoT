import math

def calculate_and_print_distances():
    """
    Calculates and prints the possible normalized distances (r <= 3) between
    hard spheres packed on a 2D hyperplane (hexagonal lattice).
    """
    print("For hard spheres on a hyperplane, the most compact arrangement is a 2D hexagonal lattice.")
    print("The normalized distance r is the distance between sphere centers divided by the sphere diameter d.\n")
    
    # Using the formula for distances in a hexagonal lattice: r^2 = i^2 + j^2 + i*j
    # where i and j are integer coordinates.
    # We find the distances r <= 3.

    # We store distances in a set to keep only unique values,
    # and a list of tuples to keep track of the calculation for printing.
    unique_distances_sq = set()
    distance_calculations = []
    
    # Iterate through possible integer coordinates i and j
    # A search limit of 4 is sufficient as 4^2 = 16 > 3^2
    limit = 4
    for i in range(-limit, limit + 1):
        for j in range(-limit, limit + 1):
            r_sq = i**2 + j**2 + i*j
            # We are interested in distances r <= 3, so r^2 <= 9
            if r_sq <= 9 and r_sq > 1e-9: # Exclude the origin (0,0)
                if r_sq not in unique_distances_sq:
                    unique_distances_sq.add(r_sq)
                    distance_calculations.append({
                        "r_sq": r_sq,
                        "i": i,
                        "j": j,
                        "r": math.sqrt(r_sq)
                    })

    # Sort distances in ascending order
    distance_calculations.sort(key=lambda x: x["r"])
    
    final_distances = []
    print("The possible distances r are:\n")

    for calc in distance_calculations:
        i, j, r_sq, r = calc["i"], calc["j"], calc["r_sq"], calc["r"]
        final_distances.append(f"{r:.2f}")

        print(f"Distance {len(final_distances)}:")
        print(f"This distance can be found from lattice point (i={i}, j={j}).")
        print(f"The squared distance r^2 is given by the equation: r^2 = i^2 + j^2 + i*j")
        # Output each number in the equation
        print(f"r^2 = {i}^2 + {j}^2 + {i}*{j} = {i**2} + {j**2} + {i*j} = {r_sq}")
        print(f"r = sqrt({r_sq:.2f}) = {r:.2f}\n")

    # Present the final set of distances as requested
    print("The set of unique distances r <= 3 is:")
    print(", ".join(final_distances))


if __name__ == "__main__":
    calculate_and_print_distances()
    # The final answer in the required format
    final_answer = [1.00, 1.73, 2.00, 2.65, 3.00]
    print(f"\n<<<{final_answer}>>>")