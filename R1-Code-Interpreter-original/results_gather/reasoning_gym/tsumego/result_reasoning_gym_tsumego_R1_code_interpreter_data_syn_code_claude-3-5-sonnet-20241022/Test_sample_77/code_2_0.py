def check_capture(board, move_x, move_y):
    # Make a copy of the board
    new_board = [row[:] for row in board]
    new_board[move_y][move_x] = 'X'  # Place black stone
    
    def has_liberties(x, y, color, visited=None):
        if visited is None:
            visited = set()
        
        if x < 0 or x >= 11 or y < 0 or y >= 11:
            return False
        
        if board[y][x] == '.':
            return True
            
        if board[y][x] != color or (x, y) in visited:
            return False
            
        visited.add((x, y))
        
        # Check all adjacent positions
        for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            new_x, new_y = x + dx, y + dy
            if has_liberties(new_x, new_y, color, visited):
                return True
        
        return False

    # Check surrounding white groups
    captures = []
    for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
        x, y = move_x + dx, move_y + dy
        if 0 <= x < 11 and 0 <= y < 11 and board[y][x] == 'O':
            if not has_liberties(x, y, 'O'):
                captures.append((x, y))
    
    return captures

# Initialize the board
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

# Test B4 move
move = (1, 7)  # B4 coordinates
captures = check_capture(board, move[0], move[1])
print(f"Move at B4 would capture stones at: {', '.join([chr(65+x)+str(y+1) for x,y in captures])}")

# Also check the current liberties of the white group
def count_liberties(x, y, color, visited=None):
    if visited is None:
        visited = set()
    
    if x < 0 or x >= 11 or y < 0 or y >= 11:
        return 0
    
    if board[y][x] != color:
        if board[y][x] == '.':
            return 1
        return 0
    
    if (x, y) in visited:
        return 0
        
    visited.add((x, y))
    liberties = 0
    
    for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
        new_x, new_y = x + dx, y + dy
        liberties += count_liberties(new_x, new_y, color, visited)
    
    return liberties

# Check liberties of key white stones
print("\nCurrent liberties:")
key_stones = [(1, 5), (2, 4), (3, 4)]  # B6, C5, D5
for x, y in key_stones:
    liberties = count_liberties(x, y, 'O')
    print(f"White stone at {chr(65+x)}{y+1} has {liberties} liberties")