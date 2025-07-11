import collections

def apply_permutation(p, subject):
    """Applies permutation p to subject tuple. p and subject are 1-indexed."""
    return tuple(subject[i-1] for i in p)

def compose_permutations(p1, p2):
    """Composes two permutations p1 o p2. p1, p2 are 1-indexed."""
    return tuple(p1[p2[i-1]-1] for i in range(1, len(p1) + 1))

def perm_add_const(p, k, n=3):
    """Adds a constant k to a permutation, modulo n. Result is 1-indexed."""
    return tuple((val - 1 + k) % n + 1 for val in p)

def get_translation_cycles(n=3):
    """Generates horizontal translation cycle permutations."""
    base = list(range(1, n + 1))
    cycles = [(tuple(base))]
    current = list(base)
    for _ in range(n - 1):
        current = current[1:] + current[:1]
        cycles.append(tuple(current))
    return [tuple(c) for c in cycles]

def get_w0(n=3):
    """Gets the reversing permutation w0."""
    return tuple(range(n, 0, -1))

def find_equivalence_classes():
    """
    Calculates the number of equivalence classes for LHTF minimal diagrams
    under translation and rotation.
    """
    # Minimal LHTF diagrams from literature
    lhtf_diagrams = {
        # Diagram A
        ((2, 3, 1), (1, 2, 3)),
        # Diagram B
        ((1, 3, 2), (2, 1, 3)),
        # Diagram C
        ((3, 2, 1), (1, 3, 2)),
    }

    n = 3
    unclassified = set(lhtf_diagrams)
    classes = []
    
    h_cycles = get_translation_cycles(n)
    w0 = get_w0(n)

    while unclassified:
        # Start a new class with an arbitrary diagram
        start_node = unclassified.pop()
        current_class = {start_node}
        queue = collections.deque([start_node])

        while queue:
            s, t = queue.popleft()

            # Generate neighbors under all symmetries
            for k in range(1, n):
                # Vertical translations
                s_v = perm_add_const(s, k, n)
                t_v = perm_add_const(t, k, n)
                neighbor_v = (s_v, t_v)
                if neighbor_v in unclassified:
                    unclassified.remove(neighbor_v)
                    current_class.add(neighbor_v)
                    queue.append(neighbor_v)

                # Horizontal translations
                c_inv = h_cycles[n-k]
                s_h = compose_permutations(h_cycles[k], compose_permutations(s, c_inv))
                t_h = compose_permutations(h_cycles[k], compose_permutations(t, c_inv))
                neighbor_h = (s_h, t_h)
                if neighbor_h in unclassified:
                    unclassified.remove(neighbor_h)
                    current_class.add(neighbor_h)
                    queue.append(neighbor_h)
            
            # 180-degree rotation
            s_r180 = compose_permutations(w0, compose_permutations(s, w0))
            t_r180 = compose_permutations(w0, compose_permutations(t, w0))
            neighbor_r180 = (s_r180, t_r180)
            if neighbor_r180 in unclassified:
                unclassified.remove(neighbor_r180)
                current_class.add(neighbor_r180)
                queue.append(neighbor_r180)

        classes.append(current_class)

    num_classes = len(classes)
    print(f"Number of minimal LHTF diagrams: {len(lhtf_diagrams)}")
    print(f"Number of equivalence classes under translation and rotation: {num_classes}")
    print(f"The number of unique diagrams is therefore: {num_classes}")


find_equivalence_classes()