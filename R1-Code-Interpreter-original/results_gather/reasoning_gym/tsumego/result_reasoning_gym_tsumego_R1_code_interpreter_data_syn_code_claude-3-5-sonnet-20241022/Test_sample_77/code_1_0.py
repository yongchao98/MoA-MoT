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

# Test B4 move
move = (1, 7)  # B4 coordinates
captures = check_capture(board, move[0], move[1])
print(f"Move at B4 would capture stones at: {', '.join([chr(65+x)+str(y+1) for x,y in captures])}")