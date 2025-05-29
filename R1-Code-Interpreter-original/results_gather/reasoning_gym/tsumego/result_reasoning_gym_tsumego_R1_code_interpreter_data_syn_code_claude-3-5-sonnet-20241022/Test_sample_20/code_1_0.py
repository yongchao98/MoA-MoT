def print_board(board):
    print("   A B C D E F G H I J K")
    for i in range(11, 0, -1):
        row = f"{i:2d} "
        for j in range(11):
            row += board[i-1][j] + " "
        print(row)

def is_valid_move(x, y):
    return 0 <= x < 11 and 0 <= y < 11

def count_group_liberties(board, x, y, color, visited=None):
    if visited is None:
        visited = set()
    
    if not is_valid_move(x, y) or (x, y) in visited:
        return set()
    
    if board[y][x] == '.':
        return {(x, y)}
    
    if board[y][x] != color:
        return set()
    
    visited.add((x, y))
    liberties = set()
    
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        new_liberties = count_group_liberties(board, x+dx, y+dy, color, visited)
        liberties.update(new_liberties)
    
    return liberties

def analyze_capture_moves():
    board = [
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', 'X', '.', '.', 'X', '.'],
        ['.', '.', '.', 'O', '.', '.', '.', 'X', 'X', 'O', 'X'],
        ['.', '.', 'O', '.', '.', '.', '.', 'X', 'O', 'O', '.'],
        ['.', '.', '.', '.', '.', '.', 'X', '.', 'X', 'O', 'X'],
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.'],
        ['.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', 'O', '.', '.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.'],
        ['.', '.', '.', '.', '.', '.', '.', '.', 'X', 'O', '.']
    ]

    # Focus on the right side where there's a cluster of stones
    potential_moves = []
    
    # Check each empty point around white stones
    for y in range(11):
        for x in range(11):
            if board[y][x] == '.':
                # Try move
                board[y][x] = 'X'
                
                # Check if this move creates a capturing opportunity
                white_groups_in_danger = False
                for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
                    nx, ny = x+dx, y+dy
                    if is_valid_move(nx, ny) and board[ny][nx] == 'O':
                        liberties = count_group_liberties(board, nx, ny, 'O')
                        if len(liberties) == 1:  # White group is in atari
                            white_groups_in_danger = True
                            # Check if the white group is substantial (more than 2 stones)
                            group_size = len([1 for i in range(11) for j in range(11) 
                                           if board[i][j] == 'O' and 
                                           (j, i) in count_group_liberties(board, nx, ny, 'O', set())])
                            if group_size >= 2:
                                col = chr(ord('A') + x)
                                row = 11 - y
                                potential_moves.append((f"{col}{row}", group_size))
                
                board[y][x] = '.'

    # Sort moves by the size of groups they threaten
    potential_moves.sort(key=lambda x: x[1], reverse=True)
    print("Potential capturing moves (move, threatened group size):", potential_moves)

analyze_capture_moves()