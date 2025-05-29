def create_board():
    board = [['.'] * 10 for _ in range(10)]
    # Add black stones (X)
    black_positions = [(9,1), (9,3), (9,4), (8,1), (8,3), (7,0), (7,8), 
                      (6,0), (5,0), (1,5), (2,8)]
    # Add white stones (O)
    white_positions = [(9,1), (9,2), (9,3), (8,2), (8,5), (6,6), (5,7), (4,8)]
    
    for x, y in black_positions:
        board[x][y] = 'X'
    for x, y in white_positions:
        board[x][y] = 'O'
    return board

def get_liberties(board, x, y, visited=None):
    if visited is None:
        visited = set()
    
    if x < 0 or x >= 10 or y < 0 or y >= 10:
        return set()
    
    if board[x][y] == '.':
        return {(x,y)}
    
    if (x,y) in visited:
        return set()
    
    visited.add((x,y))
    stone_type = board[x][y]
    liberties = set()
    
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < 10 and 0 <= new_y < 10:
            if board[new_x][new_y] == '.':
                liberties.add((new_x, new_y))
            elif board[new_x][new_y] == stone_type:
                liberties.update(get_liberties(board, new_x, new_y, visited))
    
    return liberties

def find_best_capture():
    board = create_board()
    white_groups = []
    visited = set()
    
    # Find white groups
    for i in range(10):
        for j in range(10):
            if board[i][j] == 'O' and (i,j) not in visited:
                group = set()
                stack = [(i,j)]
                while stack:
                    x, y = stack.pop()
                    if (x,y) not in visited:
                        visited.add((x,y))
                        group.add((x,y))
                        for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
                            new_x, new_y = x + dx, y + dy
                            if 0 <= new_x < 10 and 0 <= new_y < 10 and board[new_x][new_y] == 'O':
                                stack.append((new_x, new_y))
                if group:
                    white_groups.append(group)
    
    # Find liberties for each group
    best_move = None
    max_capture = 0
    
    for group in white_groups:
        liberties = set()
        for x, y in group:
            liberties.update(get_liberties(board, x, y))
        
        if len(liberties) == 1:
            liberty = liberties.pop()
            if len(group) > max_capture:
                max_capture = len(group)
                best_move = liberty

    if best_move:
        print(f"Best move: {chr(65+best_move[1])}{10-best_move[0]}")
        print(f"Captures {max_capture} stones")

find_best_capture()