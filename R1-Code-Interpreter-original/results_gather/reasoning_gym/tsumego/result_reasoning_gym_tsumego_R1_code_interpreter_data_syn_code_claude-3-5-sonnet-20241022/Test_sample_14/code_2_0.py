def analyze_move(board, move_x, move_y):
    # Make a copy of the board
    new_board = [row[:] for row in board]
    new_board[move_y][move_x] = 'X'  # Place black stone
    
    def count_liberties(x, y, color, visited=None):
        if visited is None:
            visited = set()
        
        if x < 0 or x >= 9 or y < 0 or y >= 9:
            return set()
        
        if new_board[y][x] != color:
            if new_board[y][x] == '.':
                return {(x, y)}
            return set()
        
        if (x, y) in visited:
            return set()
        
        visited.add((x, y))
        liberties = set()
        
        for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            new_x, new_y = x + dx, y + dy
            liberties.update(count_liberties(new_board, new_x, new_y, color, visited))
            if 0 <= new_x < 9 and 0 <= new_y < 9 and new_board[new_y][new_x] == '.':
                liberties.add((new_x, new_y))
        
        return liberties

    # Check white group at A3
    white_liberties = count_liberties(new_board, 0, 6, 'O')
    print(f"After move at {chr(move_x+65)}{9-move_y}, white liberties:", 
          [(chr(x+65), 9-y) for x, y in white_liberties])

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

# Test B4 move
analyze_move(board, 1, 5)  # B4 in board coordinates