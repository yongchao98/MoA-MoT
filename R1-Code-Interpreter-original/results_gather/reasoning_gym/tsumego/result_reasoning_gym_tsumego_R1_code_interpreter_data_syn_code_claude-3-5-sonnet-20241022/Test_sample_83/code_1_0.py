def create_board():
    board = [['.'] * 11 for _ in range(11)]
    # Place black stones (X)
    black_positions = [(8,2), (6,3), (6,4), (6,5), (4,5), (3,5), (3,2), (10,0)]
    # Place white stones (O)
    white_positions = [(9,2), (7,1), (6,2), (6,2), (5,4), (5,5), (5,8), (4,0)]
    
    for x, y in black_positions:
        board[x][y] = 'X'
    for x, y in white_positions:
        board[x][y] = 'O'
    return board

def find_group(board, x, y, color, visited=None):
    if visited is None:
        visited = set()
    
    if (x, y) in visited or x < 0 or y < 0 or x >= 11 or y >= 11 or board[x][y] != color:
        return set()
    
    group = {(x, y)}
    visited.add((x, y))
    
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        new_x, new_y = x + dx, y + dy
        group.update(find_group(board, new_x, new_y, color, visited))
    
    return group

def count_liberties(board, group):
    liberties = set()
    for x, y in group:
        for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
            new_x, new_y = x + dx, y + dy
            if 0 <= new_x < 11 and 0 <= new_y < 11 and board[new_x][new_y] == '.':
                liberties.add((new_x, new_y))
    return liberties

def analyze_position():
    board = create_board()
    white_groups = []
    checked = set()
    
    # Find all white groups
    for i in range(11):
        for j in range(11):
            if board[i][j] == 'O' and (i,j) not in checked:
                group = find_group(board, i, j, 'O')
                checked.update(group)
                if group:
                    liberties = count_liberties(board, group)
                    white_groups.append((group, liberties))
    
    # Analyze each white group
    best_move = None
    max_captures = 0
    
    for group, liberties in white_groups:
        if len(liberties) == 1:
            liberty = list(liberties)[0]
            if len(group) > max_captures:
                max_captures = len(group)
                best_move = liberty
    
    if best_move:
        # Convert to Go coordinates (A-K, 1-11)
        col = chr(best_move[1] + ord('A'))
        row = 11 - best_move[0]
        print(f"Best move: {col}{row} capturing {max_captures} stones")

analyze_position()