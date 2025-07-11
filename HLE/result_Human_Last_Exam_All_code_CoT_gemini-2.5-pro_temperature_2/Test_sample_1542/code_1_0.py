import collections

def count_quadratic_form_classes():
    """
    Counts the number of equivalence classes of quadratic forms in two variables
    over the ring R=Z/8Z.
    """
    N = 8
    
    # Step 1: Generate the group GL(2, Z/N)
    gl2_zN = []
    units = {i for i in range(N) if collections.gcd(i, N) == 1}
    for p in range(N):
        for q in range(N):
            for r in range(N):
                for s in range(N):
                    det = (p * s - q * r) % N
                    if det in units:
                        gl2_zN.append(((p, q), (r, s)))

    # Step 2: Generate all quadratic forms
    all_forms = []
    for a in range(N):
        for b in range(N):
            for c in range(N):
                all_forms.append((a, b, c))

    # Step 3 & 4: Count the orbits
    seen_forms = set()
    num_classes = 0
    
    for form in all_forms:
        if form in seen_forms:
            continue
        
        num_classes += 1
        
        # Compute the orbit of the current form
        orbit = set()
        queue = collections.deque([form])
        orbit.add(form)

        while queue:
            current_form = queue.popleft()
            a, b, c = current_form

            for T in gl2_zN:
                p, q = T[0]
                r, s = T[1]

                # Apply the transformation
                # a' = ap^2 + bpr + cr^2
                # b' = 2apq + b(ps+qr) + 2crs
                # c' = aq^2 + bqs + cs^2
                det = (p * s - q * r)
                
                new_a = (a*p*p + b*p*r + c*r*r) % N
                new_b = (2*a*p*q + b*det + 2*c*r*s) % N
                new_c = (a*q*q + b*q*s + c*s*s) % N
                
                new_form = (new_a, new_b, new_c)
                
                if new_form not in orbit:
                    orbit.add(new_form)
                    queue.append(new_form)

        seen_forms.update(orbit)
        
    print(f"The number of equivalence classes of quadratic forms is: {num_classes}")

count_quadratic_form_classes()