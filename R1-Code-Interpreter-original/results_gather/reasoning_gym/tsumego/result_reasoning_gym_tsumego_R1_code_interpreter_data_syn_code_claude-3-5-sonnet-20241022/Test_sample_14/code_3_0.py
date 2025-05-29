def analyze_position(board):
    def get_group_liberties(x, y, color, visited=None):
        if visited is None:
            visited = set()
        
        if x < 0 or x >= 9 or y < 0 or y >= 9:
            return set()
        
        if board[y][x] != color:
            if board[y][x] == '.':
                return {(x, y)}
            return set()
        
        if (x, y) in visited:
            return set()
        
        visited.add((x, y))
        liberties = set()
        
        for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            new_x, new_y = x + dx, y + dy
            liberties.update(get_group_liberties(new_x, new_y, color, visited))
            if 0 <= new_x < 9 and 0 <= new_y < 9 and board[new_y][new_x] == '.':
                liberties.add((new_x, new_y))
        
        return liberties

    def test_move(x, y):
        # Make a copy of the board
        test_board = [row[:] for row in board]
        if test_board[y][x] != '.':
            return False
        
        test_board[y][x] = 'X'
        
        # Check white group's liberties after the move
        liberties = get_group_liberties(0, 6, 'O')  # Starting from A3 (white stone)
        print(f"After move at {chr(x+65)}{9-y}, white liberties:", 
              [(chr(x+65), 9-y) for x, y in liberties])
        return len(liberties)

    # Initialize the board
    board = [
        ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', 'O', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.', '.', 'O'],
        ['.', '.', '.', '.', '.', '.', '.', '.', 'X'],
        ['.', '.', '.', 'X', '.', '.', '.', '.', '.'],
        ['O', '.', 'X', 'O', 'X', '.', '.', '.', 'X'],
        ['O', 'X', 'O', 'O', '.', '.', '.', '.', '.'],
        ['.', 'X', 'X', 'O', 'X', 'X', '.', '.', '.'],
        ['X', 'X', '.', 'X', '.', 'O', '.', '.', '.']
    ]

    # Test critical points
    print("Testing A5:")
    test_move(0, 4)
    print("\nTesting B6:")
    test_move(1, 3)
    print("\nTesting A8:")
    test_move(0, 1)

analyze_position(board)