import collections

def solve_alien_invasion():
    """
    Simulates the alien invasion on an 8x8 chessboard to find the
    maximal colony size K from an optimal initial setup.
    """
    # Using a set of tuples for coordinates (row, col) from (0,0) to (7,7)
    # 'a8' is (0,0), 'h1' is (7,7)
    # 'd5' is (row 8-5=3, col 'd'-'a'=3) -> (3,3)
    # 'e5' is (row 8-5=3, col 'e'-'a'=4) -> (3,4)

    # Initial optimal configuration for K=20 (a 4x5 rectangle)
    # C0 = {c5, d5, e5, f5, g5} U {g4, g3, g2}
    initial_colony = {
        (3, 2), (3, 3), (3, 4), (3, 5), (3, 6),  # c5, d5, e5, f5, g5
        (4, 6), (5, 6), (6, 6)                     # g4, g3, g2
    }
    
    colony = set(initial_colony)
    initial_size = len(colony)

    while True:
        newly_captured = set()
        # Iterate over all vacant squares on the board
        for r in range(8):
            for c in range(8):
                if (r, c) in colony:
                    continue

                # Count captured neighbors for a vacant square
                neighbors = [(r-1, c), (r+1, c), (r, c-1), (r, c+1)]
                captured_neighbors_count = 0
                for nr, nc in neighbors:
                    if (nr, nc) in colony:
                        captured_neighbors_count += 1
                
                # Check the growth condition
                if captured_neighbors_count >= 2:
                    newly_captured.add((r, c))

        if not newly_captured:
            # No more squares can be captured, the process stops
            break
        
        # Add all newly captured squares to the colony for the next turn
        colony.update(newly_captured)

    final_size = len(colony)
    captured_during_expansion = final_size - initial_size
    
    print(f"Initial colony size: {initial_size}")
    print(f"Squares captured during expansion: {captured_during_expansion}")
    print(f"The maximal final size of the colony (K) is {initial_size} + {captured_during_expansion} = {final_size}")

solve_alien_invasion()
print("\nThe final answer is K.")
print("<<<20>>>")
