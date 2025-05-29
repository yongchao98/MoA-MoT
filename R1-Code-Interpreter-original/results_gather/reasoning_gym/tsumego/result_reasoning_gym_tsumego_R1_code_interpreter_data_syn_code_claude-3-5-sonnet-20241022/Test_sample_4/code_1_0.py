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

def count_liberties(board, pos, checked=None):
    if checked is None:
        checked = set()
    
    if board[pos] == '.':
        return 1
    
    if pos in checked:
        return 0
        
    checked.add(pos)
    liberties = 0
    stone_type = board[pos]
    
    for adj in get_adjacent(pos):
        if board[adj] == '.':
            liberties += 1
        elif board[adj] == stone_type:
            liberties += count_liberties(board, adj, checked)
    
    return liberties

def find_vulnerable_stones(board):
    vulnerable = []
    for pos, stone in board.items():
        if stone == 'O':  # Check white stones
            liberties = count_liberties(board, pos)
            if liberties == 1:
                # Find the empty adjacent point
                for adj in get_adjacent(pos):
                    if board[adj] == '.':
                        vulnerable.append(adj)
    return list(set(vulnerable))

board = create_board()
capturing_moves = find_vulnerable_stones(board)
print("Possible capturing moves:", capturing_moves)

# Additionally, let's print the liberties of each white group
for pos, stone in board.items():
    if stone == 'O':
        liberties = count_liberties(board, pos)
        print(f"White stone at {pos} has {liberties} liberties")