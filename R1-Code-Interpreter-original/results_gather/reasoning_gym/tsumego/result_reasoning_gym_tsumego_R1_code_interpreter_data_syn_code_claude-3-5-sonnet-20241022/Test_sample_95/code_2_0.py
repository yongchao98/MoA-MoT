def find_captures():
    # Initialize the board (9x9)
    board = [
        [".", ".", ".", ".", ".", ".", ".", ".", "."],
        [".", ".", ".", ".", "O", ".", ".", "O", "."],
        ["X", ".", ".", ".", "X", "O", "X", ".", "X"],
        [".", ".", "X", "X", "O", "O", "X", ".", "."],
        [".", ".", ".", "X", "O", ".", "X", ".", "."],
        [".", ".", ".", ".", "X", ".", "O", ".", "."],
        [".", ".", ".", ".", "O", ".", "O", ".", "."],
        [".", ".", ".", ".", ".", ".", ".", ".", "."],
        [".", ".", ".", ".", ".", ".", ".", ".", "."]
    ]

    def is_valid_pos(x, y):
        return 0 <= x < 9 and 0 <= y < 9

    def check_capture(x, y, test_board):
        if not is_valid_pos(x, y) or test_board[x][y] != 'O':
            return 0
        
        visited = set()
        group = set()
        has_liberty = False
        
        def flood_fill(cx, cy):
            nonlocal has_liberty
            if not is_valid_pos(cx, cy) or (cx, cy) in visited:
                return
            visited.add((cx, cy))
            
            if test_board[cx][cy] == '.':
                has_liberty = True
                return
            if test_board[cx][cy] != 'O':
                return
                
            group.add((cx, cy))
            for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
                flood_fill(cx + dx, cy + dy)
        
        flood_fill(x, y)
        return len(group) if not has_liberty else 0

    captures_by_move = []
    
    # Try each empty position
    for i in range(9):
        for j in range(9):
            if board[i][j] == '.':
                # Make a copy of the board for testing
                test_board = [row[:] for row in board]
                test_board[i][j] = 'X'
                
                total_captures = 0
                # Check all adjacent positions for captures
                for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
                    nx, ny = i + dx, j + dy
                    if is_valid_pos(nx, ny):
                        captures = check_capture(nx, ny, test_board)
                        total_captures += captures
                
                if total_captures > 0:
                    move = f"{chr(ord('A') + j)}{9-i}"
                    captures_by_move.append((move, total_captures))
    
    # Sort moves by number of captures
    captures_by_move.sort(key=lambda x: x[1], reverse=True)
    
    # Print results
    if captures_by_move:
        for move, stones in captures_by_move:
            print(f"Move {move}: captures {stones} stones")
    else:
        print("No capturing moves found")

find_captures()