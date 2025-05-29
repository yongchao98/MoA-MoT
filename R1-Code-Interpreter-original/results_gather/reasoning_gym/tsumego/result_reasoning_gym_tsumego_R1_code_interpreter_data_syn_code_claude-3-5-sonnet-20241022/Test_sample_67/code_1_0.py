def create_board():
    board = {}
    board_str = """   A B C D E F G H I J K L
12 . . . . . . . . . . . .
11 . . . . O . . . . X . .
10 . . . . . . . . . . . .
 9 . . . X . O . . . X . .
 8 . . . . . . . . X O X .
 7 . O . . . . . X O O . .
 6 . . . O . . . . X O X O
 5 . . . . . X . . . X . .
 4 . . . . . . . . . . O .
 3 . . . . . . . . . . . .
 2 . . O X . . . . . . . .
 1 . . . . . . . . . . . ."""
    
    rows = board_str.strip().split('\n')[1:]
    for row_num, row in enumerate(rows):
        row_pieces = row[3:].split()
        for col_num, piece in enumerate(row_pieces):
            if piece != '.':
                board[(col_num, 12-row_num)] = piece
    return board

def find_group(board, start_x, start_y):
    color = board.get((start_x, start_y))
    if not color:
        return set()
    
    group = set([(start_x, start_y)])
    to_check = [(start_x, start_y)]
    
    while to_check:
        x, y = to_check.pop()
        for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
            new_x, new_y = x + dx, y + dy
            if (new_x, new_y) in board and board[(new_x, new_y)] == color and (new_x, new_y) not in group:
                group.add((new_x, new_y))
                to_check.append((new_x, new_y))
    return group

def get_group_liberties(board, group):
    liberties = set()
    for x, y in group:
        for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
            new_x, new_y = x + dx, y + dy
            if 0 <= new_x < 12 and 0 <= new_y < 12:
                if (new_x, new_y) not in board:
                    liberties.add((new_x, new_y))
    return liberties

def find_capturable_groups():
    board = create_board()
    checked = set()
    captures = []
    
    for (x, y), piece in board.items():
        if piece == 'O' and (x, y) not in checked:  # White stone
            group = find_group(board, x, y)
            checked.update(group)
            liberties = get_group_liberties(board, group)
            if len(liberties) == 1:
                liberty = liberties.pop()
                captures.append(liberty)
    
    cols = 'ABCDEFGHIJKL'
    for x, y in captures:
        print(f"Can capture group by playing at: {cols[x]}{y+1}")

find_capturable_groups()