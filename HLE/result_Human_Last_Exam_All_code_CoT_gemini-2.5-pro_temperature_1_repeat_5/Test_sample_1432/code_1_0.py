def solve():
    """
    Calculates the minimum number of chess pieces to remove from a 7x8 board
    to prevent any 5-in-a-row sequences.
    The solution is based on removing pieces at coordinates (r, c)
    that satisfy the condition (2*r + c) % 5 == 0, where r and c are 1-indexed.
    """
    rows = 7
    cols = 8
    
    removals_per_row = []
    total_removed = 0
    
    # Iterate through each row (1 to 7)
    for r in range(1, rows + 1):
        count_in_row = 0
        # Iterate through each column (1 to 8)
        for c in range(1, cols + 1):
            # Check if the piece should be removed
            if (2 * r + c) % 5 == 0:
                count_in_row += 1
        removals_per_row.append(count_in_row)
        total_removed += count_in_row
        
    # Print the breakdown of the final calculation
    equation_parts = [str(n) for n in removals_per_row]
    print(f"Number of pieces removed from each row: {', '.join(equation_parts)}")
    print(f"The final calculation is an addition of these numbers:")
    # Output each number in the final equation
    for i, part in enumerate(equation_parts):
        print(part, end="")
        if i < len(equation_parts) - 1:
            print(" + ", end="")
    print(f" = {total_removed}")

solve()
