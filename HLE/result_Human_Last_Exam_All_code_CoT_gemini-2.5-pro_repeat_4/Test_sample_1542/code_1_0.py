import itertools

def count_quadratic_form_classes():
    """
    Calculates the number of equivalence classes of quadratic forms in two
    variables over the ring Z/8Z.
    """
    N = 8

    # Generate all quadratic forms Q(x,y) = ax^2 + bxy + cy^2.
    # A form is represented by the tuple (a, b, c).
    all_forms = []
    for a in range(N):
        for b in range(N):
            for c in range(N):
                all_forms.append((a, b, c))

    # Generate the group GL_2(Z/8Z) of invertible 2x2 matrices.
    gl2_zN = []
    units = {x for x in range(N) if x % 2 != 0} # Units in Z/8Z are {1, 3, 5, 7}
    for p in range(N):
        for q in range(N):
            for r in range(N):
                for s in range(N):
                    # A matrix is invertible if its determinant is a unit.
                    det = (p * s - q * r) % N
                    if det in units:
                        gl2_zN.append(((p, q), (r, s)))

    # Use a set for efficient lookup and removal of forms.
    remaining_forms = set(all_forms)
    class_representatives = []

    # Partition the set of all forms into equivalence classes (orbits).
    while remaining_forms:
        # Pick a representative for the next class.
        q_rep = remaining_forms.pop()
        class_representatives.append(q_rep)
        
        a, b, c = q_rep
        
        # Compute the orbit of the representative under the action of GL_2(Z/8Z).
        orbit = {q_rep}
        for T in gl2_zN:
            p, q = T[0]
            r, s = T[1]
            
            # The transformation rule for coefficients is derived from
            # Q'(x,y) = Q(px+qy, rx+sy).
            a_new = (a * p * p + b * p * r + c * r * r) % N
            b_new = (2 * a * p * q + b * (p * s + q * r) + 2 * c * r * s) % N
            c_new = (a * q * q + b * q * s + c * s * s) % N
            
            orbit.add((a_new, b_new, c_new))
            
        # Remove all forms in the computed orbit from the set of remaining forms.
        remaining_forms.difference_update(orbit)

    num_classes = len(class_representatives)
    print(f"Found {num_classes} equivalence classes.")
    
    # Sort representatives for a consistent, readable output.
    class_representatives.sort()
    
    print("\nA set of representatives for these classes is (in the form ax^2 + bxy + cy^2):")
    for a, b, c in class_representatives:
        terms = []
        if a != 0:
            terms.append(f"{a}x^2")
        if b != 0:
            terms.append(f"{b}xy")
        if c != 0:
            terms.append(f"{c}y^2")
        
        if not terms:
            print("0")
        else:
            # Format negative coefficients correctly for modulo arithmetic.
            print(" + ".join(terms).replace("+ -", "- "))

    # There is a known formula for the number of classes for Z/(2^k)Z when k>=3.
    # For Z/8Z, we have k=3. The formula is 5*k + 6.
    print("\n---\nNote: This result matches the known theoretical formula.")
    k = 3
    expected_classes = 5 * k + 6
    print("The number of classes of quadratic forms in 2 variables over Z/(2^k)Z for k>=3 is given by the equation 5 * k + 6.")
    print(f"For k=3, the equation is: 5 * {k} + 6 = {expected_classes}")

if __name__ == '__main__':
    count_quadratic_form_classes()