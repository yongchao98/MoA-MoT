def create_board():
    board = {}
    for i in range(1, 12):
        for c in 'ABCDEFGHIJK':
            board[f"{c}{i}"] = '.'
    
    # Set white stones (O)
    white_stones = ['A11', 'H9', 'G6', 'H5', 'I5', 'D4', 'H4', 'I4', 'I3']
    for pos in white_stones:
        board[pos] = 'O'
    
    # Set black stones (X)
    black_stones = ['B7', 'I7', 'I6', 'C5', 'J5', 'G4', 'J4', 'H3', 'J3', 'I2', 'C1', 'G1']
    for pos in black_stones:
        board[pos] = 'X'
    
    return board

def get_adjacent(pos):
    col, row = pos[0], int(pos[1:])
    adj = []
    if row > 1:
        adj.append(f"{col}{row-1}")
    if row < 11:
        adj.append(f"{col}{row+1}")
    if col > 'A':
        adj.append(f"{chr(ord(col)-1)}{row}")
    if col < 'K':
        adj.append(f"{chr(ord(col)+1)}{row}")
    return adj

def check_capture_possibility(board):
    white_groups = []
    checked = set()
    
    # Find white groups
    for pos, stone in board.items():
        if stone == 'O' and pos not in checked:
            group = set()
            liberties = set()
            stack = [pos]
            
            while stack:
                current = stack.pop()
                if current not in group:
                    group.add(current)
                    checked.add(current)
                    
                    for adj in get_adjacent(current):
                        if board[adj] == 'O' and adj not in group:
                            stack.append(adj)
                        elif board[adj] == '.':
                            liberties.add(adj)
            
            white_groups.append((group, liberties))
    
    # Find potential capturing moves
    capturing_moves = []
    for group, liberties in white_groups:
        if len(liberties) == 1:
            capturing_moves.append(list(liberties)[0])
    
    return capturing_moves

board = create_board()
capturing_moves = check_capture_possibility(board)
print("Possible capturing moves:", capturing_moves)