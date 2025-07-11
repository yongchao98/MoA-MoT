import collections

def get_neighbors(r, c):
    """Returns neighbors for a cell (r, c) on an 8x8 board."""
    neighbors = []
    if r > 0: neighbors.append((r - 1, c))
    if r < 7: neighbors.append((r + 1, c))
    if c > 0: neighbors.append((r, c - 1))
    if c < 7: neighbors.append((r, c + 1))
    return neighbors

def to_coords(s):
    """Converts chess notation like 'd5' to (row, col) tuple."""
    col = ord(s[0]) - ord('a')
    row = 8 - int(s[1])
    return (row, col)

def main():
    """
    Sets up the optimal initial configuration and runs the simulation
    to find the maximal size of the alien colony.
    """
    # The chosen initial configuration consists of 4 pairs of squares,
    # including the fixed ones (d5, e5). These pairs are strategically
    # placed to seed the growth of a 6x6 square.
    # The initial perimeter is 4 * 6 = 24.
    # The target 6x6 square also has a perimeter of 2 * (6+6) = 24.
    initial_config = [
        to_coords('d5'), to_coords('e5'),  # Fixed pair
        to_coords('c4'), to_coords('f4'),  # Other 6 squares forming 3 pairs
        to_coords('c7'), to_coords('f7'),
        to_coords('d2'), to_coords('e2'),
    ]

    colony = set(initial_config)
    initial_squares = len(colony)
    
    print(f"The fixed squares are d5 and e5.")
    print(f"An optimal placement for the other 6 squares is c4, f4, c7, f7, d2, e2.")
    print(f"This initial colony of {initial_squares} squares will now expand.")
    
    turn = 0
    added_squares_count = 0
    
    while True:
        # Find all vacant squares with >= 2 captured neighbors
        eligible_squares = collections.defaultdict(list)
        for r in range(8):
            for c in range(8):
                if (r, c) not in colony:
                    neighbor_count = 0
                    for nr, nc in get_neighbors(r, c):
                        if (nr, nc) in colony:
                            neighbor_count += 1
                    if neighbor_count >= 2:
                        eligible_squares[neighbor_count].append((r, c))

        # If no more squares can be captured, the process stops
        if not eligible_squares:
            break
        
        # Optimal strategy: capture a square with the minimum number of neighbors (k)
        # to slow the perimeter decrease.
        min_k = min(eligible_squares.keys())
        
        # Sort candidates for deterministic simulation
        eligible_squares[min_k].sort()
        
        # Capture one square
        square_to_add = eligible_squares[min_k][0]
        colony.add(square_to_add)
        added_squares_count += 1
        turn += 1

    final_size = len(colony)
    print(f"\nThe expansion continued for {turn} turns.")
    print(f"The final size of the colony is {final_size} squares.")
    print("The maximal size K is determined by the sum of initial and captured squares.")
    print(f"Final equation: {final_size} = {initial_squares} + {added_squares_count}")


if __name__ == '__main__':
    main()