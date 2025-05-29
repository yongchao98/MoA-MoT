def analyze_position(board, move):
    def print_liberties_after_move(x, y):
        # Make a copy of board and place the move
        new_board = [row[:] for row in board]
        new_board[y][x] = 'X'
        
        # Check liberties of adjacent white stones
        white_stones = set()
        for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nx, ny = x + dx, y + dy
            if 0 <= nx < 11 and 0 <= ny < 11 and new_board[ny][nx] == 'O':
                white_stones.add((nx, ny))
        
        # For each white stone, find its group and count liberties
        for wx, wy in white_stones:
            group = find_connected_stones(new_board, wx, wy, 'O')
            liberties = count_group_liberties(new_board, group)
            print(f"After B4, white group at {chr(65+wx)}{wy+1} has {liberties} liberties")

    def find_connected_stones(board, x, y, color, visited=None):
        if visited is None:
            visited = set()
        
        if (x, y) in visited or x < 0 or x >= 11 or y < 0 or y >= 11 or board[y][x] != color:
            return visited
        
        visited.add((x, y))
        for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            find_connected_stones(board, x + dx, y + dy, color, visited)
        return visited

    def count_group_liberties(board, group):
        liberties = set()
        for x, y in group:
            for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nx, ny = x + dx, y + dy
                if 0 <= nx < 11 and 0 <= ny < 11 and board[ny][nx] == '.':
                    liberties.add((nx, ny))
        return len(liberties)

    print_liberties_after_move(move[0], move[1])

# Initialize board
board = [
    ['.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 11
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 10
    ['.', '.', '.', '.', '.', '.', 'O', '.', '.', '.', '.'],  # 9
    ['.', '.', '.', '.', '.', '.', '.', 'X', '.', '.', '.'],  # 8
    ['.', '.', '.', 'X', '.', '.', '.', '.', 'O', '.', '.'],  # 7
    ['.', 'O', 'X', 'O', 'X', '.', '.', '.', '.', '.', '.'],  # 6
    ['.', 'X', 'O', 'O', '.', 'O', '.', '.', '.', '.', '.'],  # 5
    ['.', '.', 'X', 'O', 'O', '.', '.', '.', 'X', '.', '.'],  # 4
    ['.', 'O', '.', 'X', '.', '.', '.', '.', '.', '.', '.'],  # 3
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.', '.', '.'],  # 2
    ['O', '.', 'O', '.', '.', '.', '.', '.', '.', '.', '.']   # 1
]

# Analyze B4 move
analyze_position(board, (1, 7))