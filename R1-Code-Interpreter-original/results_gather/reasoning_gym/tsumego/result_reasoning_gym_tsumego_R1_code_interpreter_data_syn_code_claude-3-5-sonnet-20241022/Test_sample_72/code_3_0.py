def analyze_sequence(board):
    def check_connection(x1, y1, x2, y2, stone):
        # Check if two points can be connected
        paths = [(x1, y1)]
        visited = set()
        
        while paths:
            x, y = paths.pop(0)
            if (x, y) == (x2, y2):
                return True
            
            if (x, y) in visited:
                continue
                
            visited.add((x, y))
            
            for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
                nx, ny = x + dx, y + dy
                if 0 <= nx < 9 and 0 <= ny < 9:
                    if board[ny][nx] == stone or board[ny][nx] == '.':
                        paths.append((nx, ny))
        return False

    # Test G6 move
    test_board = [row[:] for row in board]
    test_board[3][6] = 'X'  # G6
    
    # Check if white groups can still connect
    can_connect_after_g6 = check_connection(5, 5, 6, 2, 'O')  # F4 to G7
    
    print("After Black G6:")
    print("Can white groups connect?", can_connect_after_g6)
    
    # Reset and test G5
    test_board = [row[:] for row in board]
    test_board[4][6] = 'X'  # G5
    
    can_connect_after_g5 = check_connection(5, 5, 6, 2, 'O')  # F4 to G7
    
    print("\nAfter Black G5:")
    print("Can white groups connect?", can_connect_after_g5)

# Initialize board
board = [
    ['.', '.', 'O', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', 'X', '.', '.', '.', 'X', 'O', 'X'],
    ['X', 'X', '.', '.', '.', 'X', 'O', 'O', '.'],
    ['.', '.', '.', 'X', 'X', '.', 'X', 'O', 'X'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', 'O', '.', '.', 'O', '.', '.', '.'],
    ['X', '.', '.', '.', '.', 'O', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.']
]

analyze_sequence(board)