def print_board(board):
    # Print for debugging
    for i, row in enumerate(board):
        print(f"{12-i:2d} {' '.join(row)}")
    print("   A B C D E F G H I J K L")

def create_board():
    # Create empty board
    board = [['.'] * 12 for _ in range(12)]
    
    # Input the exact board state as shown
    board_state = """
   A B C D E F G H I J K L
12 . . O . . . . X . . . .
11 . . . . . . X O X . . .
10 . . . . . X O O . . . .
 9 . . . . . . X O X O . .
 8 . . . . . . . X . . O .
 7 . X O . O . . . . . . .
 6 . . . . . . . . . . . .
 5 O . . . . . . . . . . .
 4 . . . . . . . . . . . .
 3 . . . . . . . . X . . .
 2 . . . X . . . . O . . .
 1 . . X . . . . . . . . .
"""
    # Parse the board state
    lines = board_state.strip().split('\n')[1:]  # Skip the header
    for i, line in enumerate(lines):
        if i == 0:  # Skip column labels
            continue
        row = line[4:].split()  # Skip row number and get pieces
        row_idx = 12 - int(line[:3])
        for j, piece in enumerate(row):
            board[row_idx][j] = piece
    
    return board

def get_group(board, x, y, color, visited=None):
    if visited is None:
        visited = set()
    
    if x < 0 or x >= 12 or y < 0 or y >= 12:
        return set()
    
    if board[x][y] != color:
        return set()
    
    if (x, y) in visited:
        return set()
    
    visited.add((x, y))
    group = {(x, y)}
    
    for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
        new_x, new_y = x + dx, y + dy
        group.update(get_group(board, new_x, new_y, color, visited))
    
    return group

def count_liberties(board, group):
    liberties = set()
    for x, y in group:
        for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            new_x, new_y = x + dx, y + dy
            if 0 <= new_x < 12 and 0 <= new_y < 12 and board[new_x][new_y] == '.':
                liberties.add((new_x, new_y))
    return liberties

def find_capturing_moves():
    board = create_board()
    capturing_moves = []
    
    # For each empty point
    for i in range(12):
        for j in range(12):
            if board[i][j] == '.':
                # Try black move here
                board[i][j] = 'X'
                
                # Check adjacent positions for white groups
                white_groups_checked = set()
                for di, dj in [(0,1), (0,-1), (1,0), (-1,0)]:
                    ni, nj = i + di, j + dj
                    if 0 <= ni < 12 and 0 <= nj < 12 and board[ni][nj] == 'O':
                        # Get the whole white group
                        white_group = get_group(board, ni, nj, 'O')
                        group_key = tuple(sorted(white_group))
                        
                        if group_key not in white_groups_checked:
                            white_groups_checked.add(group_key)
                            liberties = count_liberties(board, white_group)
                            if not liberties:
                                move = f"{chr(ord('A') + j)}{12-i}"
                                capturing_moves.append((move, len(white_group)))
                
                # Restore the board
                board[i][j] = '.'
    
    # Sort by number of stones captured (largest first)
    capturing_moves.sort(key=lambda x: x[1], reverse=True)
    print("Possible capturing moves (move, stones captured):", capturing_moves)

find_capturing_moves()