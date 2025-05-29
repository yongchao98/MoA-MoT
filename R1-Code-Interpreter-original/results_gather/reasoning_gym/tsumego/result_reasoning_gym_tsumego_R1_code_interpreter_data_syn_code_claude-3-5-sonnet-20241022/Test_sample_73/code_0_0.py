def create_board():
    board = {}
    board_str = """
   A B C D E F G H I J K L M
13 . . . . O . . . . . . . .
12 . . . . . . . . . . . . .
11 . O . . . . . . . X . O .
10 . . . . . . . . . . . . .
 9 . . . . . . . . . . . . .
 8 . . O . . . . . . O . O .
 7 . . . . . . . . . . . . .
 6 . . . . . . . . . . X . .
 5 . . . . X . . . X . X . .
 4 . . . . . X O X O O O X .
 3 . . . . . . . . X O X . .
 2 . . . . . . . . . X . . .
 1 . . . . . . . . . . . . ."""
    
    rows = board_str.strip().split('\n')
    for i, row in enumerate(rows[1:], 1):
        row_num = 14 - i
        pieces = row[3:].split()
        for j, piece in enumerate(pieces):
            col = chr(ord('A') + j)
            if piece != '.':
                board[f"{col}{row_num}"] = piece
    return board

def get_neighbors(pos):
    col, row = pos[0], int(pos[1:])
    neighbors = []
    if col > 'A': neighbors.append(f"{chr(ord(col)-1)}{row}")
    if col < 'M': neighbors.append(f"{chr(ord(col)+1)}{row}")
    if row > 1: neighbors.append(f"{col}{row-1}")
    if row < 13: neighbors.append(f"{col}{row+1}")
    return neighbors

def count_captures(board, pos):
    if pos in board:  # Position already occupied
        return 0
    
    board[pos] = 'X'  # Try placing black stone
    captured = set()
    checked = set()
    
    # Check each neighbor for white groups that can be captured
    for neighbor in get_neighbors(pos):
        if neighbor in board and board[neighbor] == 'O' and neighbor not in checked:
            group = set()
            has_liberty = False
            stack = [neighbor]
            while stack:
                current = stack.pop()
                if current not in group:
                    group.add(current)
                    for next_pos in get_neighbors(current):
                        if next_pos not in board:  # Found liberty
                            has_liberty = True
                            break
                        elif next_pos in board and board[next_pos] == 'O' and next_pos not in group:
                            stack.append(next_pos)
            
            if not has_liberty:
                captured.update(group)
            checked.update(group)
    
    del board[pos]  # Remove test stone
    return len(captured)

def find_best_move():
    board = create_board()
    best_move = None
    max_captures = 0
    
    # Try every empty position
    for col in [chr(i) for i in range(ord('A'), ord('N'))]:
        for row in range(1, 14):
            pos = f"{col}{row}"
            if pos not in board:
                captures = count_captures(board, pos)
                if captures > max_captures:
                    max_captures = captures
                    best_move = pos
    
    print(f"Best move: {best_move}")
    print(f"Captures: {max_captures}")

find_best_move()