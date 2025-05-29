def create_board():
    board = {}
    state = """
   A B C D E F G H I J K
11 . . . . . O X X . . .
10 . . X . O X O O X . .
 9 O . . O . X O . . . .
 8 . . . . . . X . . . .
 7 . . . . . X . . . . .
 6 . X . . . . . . . . .
 5 . O . . . . . X . . .
 4 X . . . . . . . . . .
 3 . . . . . . . . . . .
 2 . . . . . . . . . . .
 1 O . . . . . . . . . ."""
    
    lines = state.strip().split('\n')
    for i, line in enumerate(lines[1:], 1):
        row = 12 - i
        cols = line[1:].split()
        for j, piece in enumerate(cols):
            if piece != '.':
                col = chr(ord('A') + j)
                board[f"{col}{row}"] = piece
    return board

def get_adjacent(pos):
    col, row = pos[0], int(pos[1:])
    adj = []
    if col > 'A':
        adj.append(f"{chr(ord(col)-1)}{row}")
    if col < 'K':
        adj.append(f"{chr(ord(col)+1)}{row}")
    if row > 1:
        adj.append(f"{col}{row-1}")
    if row < 11:
        adj.append(f"{col}{row+1}")
    return adj

def analyze_position():
    board = create_board()
    white_groups = []
    checked = set()
    
    # Find white groups
    for pos in board:
        if board[pos] == 'O' and pos not in checked:
            group = set()
            liberties = set()
            stack = [pos]
            while stack:
                current = stack.pop()
                if current not in checked:
                    checked.add(current)
                    if board.get(current) == 'O':
                        group.add(current)
                        for adj in get_adjacent(current):
                            if adj not in checked:
                                if adj not in board:
                                    liberties.add(adj)
                                elif board.get(adj) == 'O':
                                    stack.append(adj)
            if group:
                white_groups.append((group, liberties))
    
    # Analyze each group
    for group, liberties in white_groups:
        print(f"White group at {sorted(group)} has liberties at {sorted(liberties)}")

analyze_position()