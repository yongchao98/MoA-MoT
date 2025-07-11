def solve_reversal_moves():
    """
    Calculates the minimum number of moves to reverse a sequence of 100 elements
    with the given swap operations.
    """
    N = 100  # Total number of elements
    K = 5    # The non-adjacent swap distance

    # There are K congruence classes modulo K.
    # The number of elements in each class is N / K.
    num_per_class = N // K

    print(f"The sequence has {N} elements.")
    print(f"Non-adjacent swaps of distance {K} are free.")
    print(f"This partitions the {N} positions into {K} classes of {num_per_class} elements each.")
    print("-" * 20)

    total_moves = 0
    class_costs = {}

    # Determine the moves required for each class
    # Classes are numbered 0, 1, 2, 3, 4
    for c in range(K):
        # An element at position p (1-indexed) is in class (p-1) % K.
        # Its target position is N+1-p.
        # Its target class is ((N+1-p)-1) % K = (N-p) % K.
        # For p in class c, p = c+1 (mod K).
        # Target class = (N - (c+1)) % K = (100 - c - 1) % K = (99-c) % K = (4-c) % 5
        # We use 0-indexed positions (p=0..99) for simplicity in calculation
        # Initial class: p % K -> c
        # Target pos: (N-1)-p. Target class: ((N-1)-p) % K
        # If p is in class c, p = m*K + c.
        # Target class = (N-1 - (m*K+c)) % K = (N-1-c) % K = (99-c)%5
        target_c = (N - 1 - c) % K
        
        # Distance is the minimum number of adjacent swaps needed to move between classes
        dist = min(abs(c - target_c), K - abs(c - target_c))
        
        cost = num_per_class * dist
        class_costs[c] = {'target': target_c, 'dist': dist, 'cost': cost}
        
    # We can group the calculations by swapping pairs of classes.
    # e.g., class 0 moves to class 4, and class 4 moves to class 0.
    
    equation_parts = []
    
    # Use a set to track processed classes to avoid double counting
    processed_classes = set()

    for c in range(K):
        if c in processed_classes:
            continue
        
        data = class_costs[c]
        target_c = data['target']
        dist = data['dist']

        if c == target_c:
            # This class maps to itself
            cost = data['cost']
            print(f"Moves for elements in class {c} to be rearranged within their target class {target_c}:")
            print(f"{num_per_class} elements * {dist} moves = {cost}")
            total_moves += cost
            equation_parts.append(str(cost))
            processed_classes.add(c)
        else:
            # These two classes swap places
            # Cost for c -> target_c is data['cost']
            # Cost for target_c -> c is class_costs[target_c]['cost'], which is the same.
            pair_cost = data['cost'] + class_costs[target_c]['cost']
            print(f"Moves for classes {c} and {target_c} to swap places:")
            print(f"({num_per_class} elements * {dist} moves) + ({num_per_class} elements * {dist} moves) = {pair_cost}")
            total_moves += pair_cost
            equation_parts.append(str(pair_cost))
            processed_classes.add(c)
            processed_classes.add(target_c)

    print("-" * 20)
    print(f"Final Equation: Total Moves = {' + '.join(equation_parts)}")
    print(f"The minimum number of moves is: {total_moves}")


solve_reversal_moves()
<<<120>>>