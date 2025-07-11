def solve():
    """
    Calculates the minimum number of moves to reverse a sequence of 100 elements
    with the given swap operations.
    """
    N = 100  # Total number of elements
    D = 5    # Non-adjacent swap distance

    num_classes = D
    class_size = N // num_classes

    # `dest_map[k]` will be the destination class for elements starting in class k.
    # Classes are 1-indexed (1 to 5).
    dest_map = {}
    for k in range(1, num_classes + 1):
        # The destination position for an element at `i` is `N+1-i`.
        # The destination class for an element in class `k` (where i % D == k % D)
        # is `(N+1-i) % D`. Since N=100 is a multiple of D=5, (100+1-i)%5 simplifies
        # to (1-i)%5. If i is in class k, i=k (mod 5). So dest class is (1-k) mod 5.
        dest = (1 - k) % num_classes
        if dest == 0:
            dest = num_classes  # Convert 0 to 5 for 1-based indexing
        dest_map[k] = dest

    total_moves = 0
    boundary_moves = []

    # Calculate the number of elements crossing each boundary k -> k+1
    for k in range(1, num_classes):  # For boundaries 1|2, 2|3, 3|4, 4|5
        crossings = 0
        # Check all source classes to the left of or at the boundary (1 to k)
        for src_class in range(1, k + 1):
            # If their destination is to the right of the boundary (k+1 to 5)
            if dest_map[src_class] > k:
                crossings += class_size
        boundary_moves.append(crossings)
        total_moves += crossings

    # Output the result with the equation
    equation_parts = [str(m) for m in boundary_moves]
    print(f"The minimum number of moves is calculated by summing the required swaps at each class boundary.")
    for i, moves in enumerate(boundary_moves):
        print(f"Moves across boundary C{i+1}|C{i+2}: {moves}")
    
    print(f"\nTotal moves = {' + '.join(equation_parts)} = {total_moves}")

solve()
<<<120>>>