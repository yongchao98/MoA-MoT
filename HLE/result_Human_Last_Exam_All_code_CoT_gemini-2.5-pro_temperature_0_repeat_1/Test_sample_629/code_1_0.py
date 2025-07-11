import collections

def get_canonical(X, O):
    """
    Finds the canonical representation of a diagram under translation.
    The canonical form is the lexicographically smallest of all n*n translations.
    A translation is a cyclic shift of columns by p and rows by q.
    """
    n = len(X)
    translations = []
    for p in range(n):  # Iterate through all column shifts
        for q in range(n):  # Iterate through all row shifts
            X_new = [0] * n
            O_new = [0] * n
            for i in range(n):
                # Apply column shift p
                col_shifted_i = (i + p) % n
                # Apply row shift q
                X_new[col_shifted_i] = (X[i] + q) % n
                O_new[col_shifted_i] = (O[i] + q) % n
            translations.append((tuple(X_new), tuple(O_new)))
    
    # Return the lexicographically smallest translation
    return min(translations)

def rotate_90_clockwise(X, O):
    """
    Rotates a grid diagram 90 degrees clockwise.
    A point at (row, col) moves to (col, n-1-row).
    This translates to a transformation on the permutation pair.
    """
    n = len(X)
    X_new = [0] * n
    O_new = [0] * n
    
    # For X
    for c, r in enumerate(X):
        new_r, new_c = c, n - 1 - r
        X_new[new_r] = new_c
        
    # For O
    for c, r in enumerate(O):
        new_r, new_c = c, n - 1 - r
        O_new[new_r] = new_c
        
    return (tuple(X_new), tuple(O_new))

def solve():
    """
    Calculates the number of unique minimal grid diagrams for the left-hand trefoil
    up to translation and rotation.
    """
    # A minimal grid diagram (n=3) for the left-hand trefoil knot.
    # Using 0-based indexing for permutations.
    initial_X = (0, 2, 1)
    initial_O = (1, 0, 2)
    
    # Use a queue for BFS and a set to store unique canonical forms
    q = collections.deque([(initial_X, initial_O)])
    initial_canonical = get_canonical(initial_X, initial_O)
    found_canonical_forms = {initial_canonical}
    
    print("Starting with the initial diagram for the left-hand trefoil:")
    print(f"  X = {list(initial_X)}, O = {list(initial_O)}")
    print(f"Its canonical form is: {initial_canonical}\n")
    
    processed_count = 0
    while q:
        current_X, current_O = q.popleft()
        processed_count += 1
        
        diagram_to_rotate = (current_X, current_O)
        
        # Generate 90, 180, 270 degree rotations
        for i in range(3):
            diagram_to_rotate = rotate_90_clockwise(diagram_to_rotate[0], diagram_to_rotate[1])
            
            # Get the canonical form of the rotated diagram
            canonical_form = get_canonical(diagram_to_rotate[0], diagram_to_rotate[1])
            
            # If this canonical form has not been seen before, add it
            if canonical_form not in found_canonical_forms:
                found_canonical_forms.add(canonical_form)
                q.append(diagram_to_rotate)

    print("The unique canonical forms found are:")
    for i, (X, O) in enumerate(sorted(list(found_canonical_forms))):
        print(f"{i+1}: X={list(X)}, O={list(O)}")

    print("\nThe total number of unique grid diagrams is the number of canonical forms found.")
    print(f"Final count: {len(found_canonical_forms)}")

solve()