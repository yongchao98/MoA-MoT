def solve_quadratic_form_classes():
    """
    Calculates and prints the equivalence classes of quadratic forms in two variables
    over the ring R=Z/8Z.
    """
    # Step 1: Generate the group GL(2, Z/8Z)
    matrices = []
    units = {1, 3, 5, 7}
    R = 8
    for p in range(R):
        for q in range(R):
            for r in range(R):
                for s in range(R):
                    det = (p * s - q * r) % R
                    if det in units:
                        matrices.append(((p, q), (r, s)))
    
    gl2_z8 = matrices

    # Step 2: Define the transformation on a form's coefficients (a,b,c)
    def transform(form, matrix):
        a, b, c = form
        p, q = matrix[0]
        r, s = matrix[1]
        
        a_new = (a * p * p + b * p * r + c * r * r) % R
        b_new = (2 * a * p * q + b * (p * s + q * r) + 2 * c * r * s) % R
        c_new = (a * q * q + b * q * s + c * s * s) % R
        
        return (a_new, b_new, c_new)

    # Step 3: Classify all 512 forms using the orbit counting algorithm
    all_forms = set()
    for a in range(R):
        for b in range(R):
            for c in range(R):
                all_forms.add((a, b, c))

    representatives = []
    
    while all_forms:
        # Pick an unclassified form as a new representative
        rep_form = next(iter(all_forms))
        representatives.append(rep_form)
        
        # Compute the orbit of this representative
        orbit = set()
        for matrix in gl2_z8:
            new_form = transform(rep_form, matrix)
            orbit.add(new_form)
            
        # Remove the entire orbit from the set of unclassified forms
        all_forms -= orbit

    # Step 4: Print the results in a structured table
    num_classes = len(representatives)
    
    class_details = []
    for rep in sorted(representatives):
        a, b, c = rep
        discriminant = (b*b - 4*a*c) % R
        
        terms = []
        if a != 0: terms.append(f"{a}x^2")
        if b != 0: terms.append(f"{b}xy")
        if c != 0: terms.append(f"{c}y^2")
        
        if not terms:
             form_str = "0"
        else:
             # This creates the equation string with each number
             form_str = " + ".join(terms)

        class_details.append({
            "rep": rep,
            "form": form_str,
            "D": discriminant
        })
    
    # Sort by discriminant then by the tuple representation for a canonical ordering
    class_details.sort(key=lambda x: (x['D'], x['rep']))
    
    print("Found {} equivalence classes. The representatives are:".format(num_classes))
    print("-" * 75)
    print("{:<8} | {:<25} | {:<25} | {}".format("Class #", "Representative (a,b,c)", "Equation", "Discriminant D=b^2-4ac"))
    print("-" * 75)

    for i, detail in enumerate(class_details):
        print("{:<8} | {:<25} | {:<25} | {}".format(
            i + 1,
            str(detail['rep']),
            detail['form'],
            detail['D']
        ))
    print("-" * 75)
    print("Total number of equivalence classes: {}".format(num_classes))


solve_quadratic_form_classes()
<<<24>>>