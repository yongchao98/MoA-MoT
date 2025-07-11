import sys

def solve():
    """
    Calculates the number of equivalence classes of quadratic forms in two
    variables over the ring R=Z/8Z.
    """
    R = 8

    # Generate the group of invertible 2x2 matrices over Z/8Z
    gl2_z8 = []
    units = {1, 3, 5, 7}
    for p in range(R):
        for q in range(R):
            for r in range(R):
                for s in range(R):
                    det = (p * s - q * r) % R
                    if det in units:
                        gl2_z8.append(((p, q), (r, s)))

    # Transformation function for a quadratic form (a,b,c) by a matrix T
    def transform(form, T):
        a, b, c = form
        p, q = T[0]
        r, s = T[1]

        a_prime = (a * p * p + b * p * r + c * r * r) % R
        b_prime = (2 * a * p * q + b * (p * s + q * r) + 2 * c * r * s) % R
        c_prime = (a * q * q + b * q * s + c * s * s) % R
        return (a_prime, b_prime, c_prime)

    # Generate all possible forms
    all_forms = []
    for a in range(R):
        for b in range(R):
            for c in range(R):
                all_forms.append((a, b, c))

    classified_forms = set()
    representatives = []
    orbit_sizes = []

    # Partition the set of all forms into equivalence classes (orbits)
    for form in all_forms:
        if form in classified_forms:
            continue

        representatives.append(form)
        
        # Compute the orbit of the current form
        orbit = set()
        for T in gl2_z8:
            new_form = transform(form, T)
            orbit.add(new_form)
        
        orbit_sizes.append(len(orbit))
        classified_forms.update(orbit)

    # Output the results
    # The "equation" is the sum of the sizes of the orbits
    orbit_sizes.sort()
    equation = " + ".join(map(str, orbit_sizes))
    total_forms = sum(orbit_sizes)
    
    print(f"The set of 512 quadratic forms is partitioned into {len(representatives)} equivalence classes.")
    print("The sizes of these classes (orbits) are:")
    print(", ".join(map(str, orbit_sizes)))
    print("\nThe sum of the orbit sizes gives the total number of forms:")
    print(f"{equation} = {total_forms}")
    print("\nThus, the number of equivalence classes of quadratic forms is:")
    print(len(representatives))

solve()
<<<22>>>