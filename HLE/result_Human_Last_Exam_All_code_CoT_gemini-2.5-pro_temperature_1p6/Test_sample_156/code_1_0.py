def solve_chess_aliens():
    """
    Solves the chess aliens problem by finding the optimal initial configuration
    and simulating the growth to find the maximal colony size K.
    """
    # Board is 8x8. We use 0-indexed (row, col) coordinates.
    # 'a'.. 'h' -> 0..7, '1'..'8' -> 0..7.
    # d5 is (4, 3)
    # e5 is (4, 4)
    
    # Based on the plan, the optimal strategy aims to fill a 5x5 rectangle,
    # leading to a maximal area of K=25.
    
    # We define an initial set S_0 of 8 squares that includes d5 and e5 and
    # forms an "L" shape for a 5x5 area (rows 1-5, cols b-f).
    # S_0 = {b5, c5, d5, e5} and {f1, f2, f3, f4}.
    initial_squares = {
        (4, 1), (4, 2), (4, 3), (4, 4),  # b5, c5, d5, e5
        (0, 5), (1, 5), (2, 5), (3, 5)   # f1, f2, f3, f4
    }
    
    captured = set(initial_squares)
    board_size = 8

    # Simulate the expansion process
    while True:
        newly_captured = set()
        # Iterate over all vacant squares
        for r in range(board_size):
            for c in range(board_size):
                if (r, c) in captured:
                    continue
                
                # Count captured neighbors for the vacant square (r,c)
                neighbor_count = 0
                for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                    nr, nc = r + dr, c + dc
                    if 0 <= nr < board_size and 0 <= nc < board_size:
                        if (nr, nc) in captured:
                            neighbor_count += 1
                
                if neighbor_count >= 2:
                    newly_captured.add((r, c))

        if not newly_captured:
            # No new squares were captured, the colony has stabilized.
            break
        
        captured.update(newly_captured)
        
    K = len(captured)
    
    print("The maximal size of the aliens' colony (K) is determined by finding the optimal initial arrangement of 8 squares (including d5 and e5) that leads to the largest possible filled area.")
    print("The strategy is to form an 'L' shape with the initial squares along the border of a target M x N rectangle.")
    print("The number of squares in such an 'L' is M + N - 2. With 8 initial squares, we have M + N - 2 = 8, so M + N = 10.")
    print("The area M*N is maximized when M=5 and N=5, giving a target area of 25.")
    print("\nAn initial set that achieves this is {b5, c5, d5, e5, f1, f2, f3, f4}.")
    print("This set will trigger a chain reaction that fills a 5x5 rectangular area.")
    
    initial_size = len(initial_squares)
    newly_captured_size = K - initial_size
    
    print("\nThe calculation for the maximal size K is:")
    print(f"Final Size (K) = Initial Size + Newly Captured Squares")
    print(f"{K} = {initial_size} + {newly_captured_size}")

solve_chess_aliens()
<<<25>>>