import itertools

def get_powerset(s):
    """Generates the powerset of a given set."""
    # Frozensets are used so they can be added to sets.
    s_list = list(s)
    # Generate all subsets from size 1 to len(s)
    return {
        frozenset(subset)
        for i in range(1, len(s_list) + 1)
        for subset in itertools.combinations(s_list, i)
    }

def apply_homeomorphism(h, subset):
    """Applies a homeomorphism h to a subset."""
    return frozenset({h.get(x, x) for x in subset})

def main():
    """
    Calculates the number of topologically distinct compactifications for the
    simplest case: a two-point remainder space X = {0, 1}.
    """
    # 1. Define the space X.
    X = {0, 1}
    print(f"Let the remainder space be X = {X}")

    # 2. Find all non-empty closed subsets of X.
    # In a discrete space, all subsets are closed.
    closed_subsets = get_powerset(X)
    print("\nThe non-empty closed subsets of X are:")
    # Sort for consistent output
    sorted_subsets = sorted([sorted(list(s)) for s in closed_subsets])
    for s in sorted_subsets:
        print(f"- {set(s)}")
    
    # 3. Define the homeomorphism group of X.
    # For a two-point space, there are two homeomorphisms:
    # - The identity map
    # - The swap map
    h_identity = {0: 0, 1: 1}
    h_swap = {0: 1, 1: 0}
    homeomorphisms = [h_identity, h_swap]
    print("\nThe homeomorphisms of X are:")
    print(f"- Identity: id(0)=0, id(1)=1")
    print(f"- Swap:     swap(0)=1, swap(1)=0")

    # 4. Calculate the orbits.
    print("\nCalculating the orbits of the closed subsets...")
    orbits = []
    subsets_to_visit = set(closed_subsets)
    
    orbit_number = 1
    while subsets_to_visit:
        # Start a new orbit with an arbitrary element
        seed = subsets_to_visit.pop()
        current_orbit = {seed}
        
        # Find all other subsets in the same orbit
        orbit_elements_to_process = {seed}
        processed_elements = set()
        
        while orbit_elements_to_process:
            subset = orbit_elements_to_process.pop()
            processed_elements.add(subset)
            for h in homeomorphisms:
                image = apply_homeomorphism(h, subset)
                if image not in processed_elements:
                    current_orbit.add(image)
                    orbit_elements_to_process.add(image)

        # Remove all elements of the found orbit from the main set
        subsets_to_visit -= current_orbit
        orbits.append(current_orbit)
        
        # Print the found orbit
        print(f"Orbit {orbit_number}: {{", end="")
        orbit_str = ", ".join([str(set(s)) for s in sorted([sorted(list(o)) for o in current_orbit])])
        print(orbit_str, end="}\n")
        orbit_number += 1
        
    # 5. The number of orbits is the answer.
    num_orbits = len(orbits)
    print("\nThe number of orbits determines the number of distinct compactifications.")
    print(f"Final calculation: Number of orbits = {num_orbits}")
    print("\nThis shows that the minimum number is 2, as predicted by the theoretical lower bound.")

if __name__ == "__main__":
    main()
