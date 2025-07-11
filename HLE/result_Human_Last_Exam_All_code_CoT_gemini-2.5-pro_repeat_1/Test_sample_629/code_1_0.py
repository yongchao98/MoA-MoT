import itertools

def perm_mult(p1, p2):
    """Computes the product of two permutations, p1 * p2."""
    n = len(p1)
    result = [0] * n
    for i in range(n):
        result[i] = p1[p2[i]]
    return tuple(result)

def find_orbits(elements, operator):
    """Partitions a set of elements into orbits under a given operator."""
    unvisited = set(elements)
    orbits = []
    while unvisited:
        start_node = unvisited.pop()
        current_orbit = {start_node}
        queue = [start_node]
        head = 0
        while head < len(queue):
            current_node = queue[head]
            head += 1
            
            # Apply the operator to find the next element in the orbit
            next_node = operator(current_node)
            if next_node not in current_orbit:
                current_orbit.add(next_node)
                queue.append(next_node)
                if next_node in unvisited:
                    unvisited.remove(next_node)
        orbits.append(sorted(list(current_orbit)))
    return orbits

def main():
    """Main function to solve the problem."""
    n = 3
    # Get all permutations of (0, 1, 2)
    all_perms = list(itertools.permutations(range(n)))

    # The characteristic permutation for the left-hand trefoil (a 3-cycle).
    p_lht = (2, 0, 1)

    # Step 1: Generate all 6 minimal grid diagrams for the left-hand trefoil.
    lht_diagrams = []
    for tau in all_perms:
        sigma = perm_mult(p_lht, tau)
        lht_diagrams.append((sigma, tau))

    print(f"Found {len(lht_diagrams)} minimal grid diagrams for the left-hand trefoil.")

    # Step 2: Define the translation operator (cyclic row shift).
    # The cycle z = (0 1 2) is represented as the permutation (1, 2, 0).
    z = (1, 2, 0)
    
    def translation_operator(diagram):
        sigma, tau = diagram
        new_sigma = perm_mult(z, sigma)
        new_tau = perm_mult(z, tau)
        return (new_sigma, new_tau)

    # Step 3: Find the orbits under the translation operator.
    orbits = find_orbits(lht_diagrams, translation_operator)
    num_orbits = len(orbits)

    print("\nThese diagrams are partitioned into orbits by 'translation' (cyclic shifts):")
    for i, orbit in enumerate(orbits):
        print(f"  Orbit {i+1} contains {len(orbit)} diagrams.")

    print("\nRotation by 90 degrees maps a left-handed trefoil to a right-handed one.")
    print("Therefore, rotation cannot reduce the number of distinct orbits of left-handed trefoils.")
    
    final_answer = num_orbits
    print(f"\nThe number of unique grid diagrams up to translation and rotation is {final_answer}.")

if __name__ == "__main__":
    main()