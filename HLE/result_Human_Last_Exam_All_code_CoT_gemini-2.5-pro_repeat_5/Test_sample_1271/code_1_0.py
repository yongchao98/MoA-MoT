import itertools

def get_homeomorphism_class(vector_set):
    """
    Determines the homeomorphism class of a set of rays from the origin.
    The set is represented by its direction vectors.
    """
    num_vectors = len(vector_set)

    if num_vectors == 0:
        # The intersection is just the origin.
        return "A single point"
    
    if num_vectors == 1:
        # The intersection is a single ray starting from the origin.
        # This is homeomorphic to [0, infinity).
        return "A ray"

    if num_vectors == 2:
        # The intersection has two rays.
        # We model vectors and their opposites abstractly with strings.
        # e.g., 'u' and '-u'.
        v1, v2 = list(vector_set)
        
        # Check if one vector string is the negative of the other.
        is_opposite = (v1.startswith('-') and v1[1:] == v2) or \
                      (v2.startswith('-') and v2[1:] == v1)

        if is_opposite:
            # The two rays are opposite, forming a line.
            # This is homeomorphic to R.
            return "A line"
        else:
            # The two rays are not opposite, forming a 'Y' shape.
            return "A bent line (Y-shape)"
            
    # The number of rays in an intersection cannot exceed 2.
    return "Invalid intersection"

def solve():
    """
    Systematically finds all possible homeomorphism classes for the
    intersection of two geodesics.
    """
    # We use strings to represent abstract, normalized, linearly independent vectors.
    # We can always find such functions in C[0,1] (e.g., normalized polynomials).
    # u, v, w, z are assumed to be linearly independent.
    # -u is the opposite vector to u.

    # Geodesics are defined by their set of 1 or 2 direction vectors.
    # Type 1: A line, with opposite direction vectors.
    line_1 = frozenset(['u', '-u'])
    line_2 = frozenset(['v', '-v'])

    # Type 2: A bent line, with two linearly independent direction vectors.
    bent_1 = frozenset(['u', 'v'])
    bent_2 = frozenset(['u', 'w'])
    bent_3 = frozenset(['w', 'z'])
    
    # We generate all possible pairs of geodesic types and find their intersections.
    geodesics_to_test = [
        line_1,    # A line L_u
        line_2,    # A different line L_v
        bent_1,    # A bent line B_{u,v}
        bent_2,    # A bent line B_{u,w} sharing one ray with B_{u,v}
        bent_3     # A bent line B_{w,z} with no rays in common with the others
    ]

    found_classes = set()

    # Iterate through all pairs of defined geodesics (including intersecting a geodesic with itself)
    for g1, g2 in itertools.product(geodesics_to_test, repeat=2):
        # The direction vectors of the intersection are the set intersection of direction vectors.
        intersection_vectors = g1.intersection(g2)
        
        # Determine the homeomorphism class of the resulting shape.
        h_class = get_homeomorphism_class(intersection_vectors)
        found_classes.add(h_class)

    print("The possible homeomorphism classes for the intersections of two geodesics are:")
    # Sort for consistent ordering
    sorted_classes = sorted(list(found_classes))
    for i, hc in enumerate(sorted_classes):
        print(f"  {i + 1}. {hc}")
    
    final_count = len(found_classes)
    print(f"\nIn total, there are {final_count} distinct homeomorphism classes.")

solve()