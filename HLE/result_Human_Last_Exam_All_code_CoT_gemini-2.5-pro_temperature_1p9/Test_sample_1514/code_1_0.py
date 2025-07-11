import collections

def find_smallest_number_of_compactifications():
    """
    This function demonstrates that the smallest number of topologically
    distinct compactifications of the ray with remainder X is 1.
    
    The number of compactifications is the number of orbits of the set of
    non-empty subcontinua of X under the action of the group of homeomorphisms of X.
    
    We choose X to be a two-point discrete space, X = {0, 1}. This space is:
    - non-degenerate (2 points)
    - locally-connected (discrete topology)
    - compact (finite)
    - metric (e.g., discrete metric d(0,1)=1)
    """

    print("Let X be the two-point space {0, 1} with the discrete topology.")
    print("This space satisfies all the required properties.")
    print("-" * 50)
    
    # In a discrete space, the only non-empty connected subsets (subcontinua)
    # are the singletons. We use frozensets because they are hashable.
    subcontinua = {frozenset([0]), frozenset([1])}
    print("The set of non-empty subcontinua of X is:", [set(s) for s in subcontinua])
    
    # For a discrete space, any bijection is a homeomorphism.
    # We can represent these as dictionaries mapping points to points.
    homeomorphisms = [
        {0: 0, 1: 1},  # Identity map
        {0: 1, 1: 0}   # Swap map
    ]
    print("The group of homeomorphisms of X consists of:")
    print("- The identity map: {0->0, 1->1}")
    print("- The swap map:    {0->1, 1->0}")
    print("-" * 50)

    # The action of a homeomorphism h on a subcontinuum s is the set {h(x) for x in s}
    def apply_homeo(h, s):
        return frozenset({h[x] for x in s})

    # We compute the orbits of the subcontinua under the action of the homeomorphisms.
    # An orbit of an element s is the set of all elements it can be mapped to.
    
    remaining = set(subcontinua)
    orbits = []
    
    while remaining:
        # Pick a representative for a new orbit
        representative = remaining.pop()
        
        # Compute the full orbit for this representative
        orbit = {representative}
        for h in homeomorphisms:
            image = apply_homeo(h, representative)
            orbit.add(image)
            
        orbits.append(orbit)
        # Remove all elements of the found orbit from the remaining set
        remaining.difference_update(orbit)
        
    print(f"We partition the set of subcontinua into orbits:")
    for i, orbit in enumerate(orbits):
        # Convert frozensets to sets for clearer printing
        printable_orbit = [set(s) for s in orbit]
        print(f"Orbit {i+1}: {printable_orbit}")
    
    print("-" * 50)
    
    num_orbits = len(orbits)
    print("The number of topologically distinct compactifications is the number of orbits.")
    print(f"The number of orbits found is: {num_orbits}")
    print("-" * 50)
    print("Therefore, the smallest number of topologically distinct compactifications is 1.")

find_smallest_number_of_compactifications()