def calculate_orbits(points, automorphisms):
    """
    Calculates the orbits of a set of points under a group of automorphisms.

    Args:
        points (list): A list of the points in the space, e.g., [0, 1, 2].
        automorphisms (list): A list of dictionaries, where each dictionary
                              represents a permutation (a homeomorphism).

    Returns:
        list: A list of sets, where each set is an orbit.
    """
    unvisited = set(points)
    orbits = []
    
    while unvisited:
        # Start a new orbit with an unvisited point
        p = unvisited.pop()
        current_orbit = {p}
        queue = [p] # Points to visit within this orbit
        head = 0
        
        while head < len(queue):
            current_point = queue[head]
            head += 1
            
            # Apply all automorphisms to find all points reachable from current_point
            for h in automorphisms:
                image = h[current_point]
                if image not in current_orbit:
                    current_orbit.add(image)
                    queue.append(image)
                    if image in unvisited:
                        unvisited.remove(image)
        
        orbits.append(current_orbit)
        
    return orbits

def main():
    """
    Main function to demonstrate the calculation for two different spaces.
    """
    print("This script calculates the number of distinct compactifications of a ray")
    print("for a given remainder space X. This number is equal to the number of orbits")
    print("of X under its group of homeomorphisms, Aut(X).\n")
    
    # --- Case 1: A homogeneous space X1 ---
    # This space is two discrete points, {0, 1}. It is nondegenerate,
    # locally connected, compact, metric, and homogeneous.
    # It is a model for the 0-sphere S^0.
    points_x1 = [0, 1]
    # The homeomorphism group consists of the identity and the swap map.
    automorphisms_x1 = [
        {0: 0, 1: 1},  # Identity
        {0: 1, 1: 0}   # Swap
    ]
    
    orbits_x1 = calculate_orbits(points_x1, automorphisms_x1)
    num_orbits_x1 = len(orbits_x1)
    
    print("--- Case 1: A Homogeneous Space (e.g., two points) ---")
    print(f"Space X1 has points: {points_x1}")
    print("Orbits are:")
    for orbit in orbits_x1:
        print(f"  {orbit}")
    print(f"The number of orbits is N(X1) = {num_orbits_x1}")
    print("This corresponds to 1 distinct compactification.")
    print("-" * 30 + "\n")

    # --- Case 2: A non-homogeneous space X2 ---
    # This space represents the interval [0, 1], modeled by a path on 3 vertices {0, 1, 2}.
    # The endpoints are 0 and 2, and the center is 1.
    # It is nondegenerate, locally connected, compact, and metric.
    points_x2 = [0, 1, 2]
    # The homeomorphisms can only swap the endpoints, but must fix the center point.
    automorphisms_x2 = [
        {0: 0, 1: 1, 2: 2},  # Identity
        {0: 2, 1: 1, 2: 0}   # Flip endpoints
    ]

    orbits_x2 = calculate_orbits(points_x2, automorphisms_x2)
    num_orbits_x2 = len(orbits_x2)

    print("--- Case 2: A Non-Homogeneous Space (e.g., an interval) ---")
    print(f"Space X2 has points: {points_x2}")
    print("Orbits are:")
    for orbit in orbits_x2:
        print(f"  {orbit}")
    print(f"The number of orbits is N(X2) = {num_orbits_x2}")
    print("This corresponds to 2 distinct compactifications.")
    print("-" * 30 + "\n")

    # --- Conclusion ---
    smallest_number = min(num_orbits_x1, num_orbits_x2)
    print("The number of compactifications depends on the choice of X.")
    print(f"Based on our analysis, the smallest number must be at least 1, and we have found a case where the number is 1.")
    print(f"Therefore, the smallest number of topologically distinct compactifications is {smallest_number}.")


if __name__ == "__main__":
    main()
<<<1>>>