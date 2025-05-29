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
 1 . . . . . . . . . . ."""
    
    rows = board_str.strip().split('\n')[1:]  # Skip the column labels
    for row_num, row in enumerate(rows):
        row_pieces = row[3:].split()  # Skip the row number
        for col_num, piece in enumerate(row_pieces):
            if piece != '.':
                board[(col_num, 12-row_num)] = piece
    return board

def get_liberties(board, x, y):
    liberties = []
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < 12 and 0 <= new_y < 12:
            if (new_x, new_y) not in board:
                liberties.append((new_x, new_y))
    return liberties

def find_captures():
    board = create_board()
    potential_captures = []
    
    # Look for white stones with few liberties
    for (x, y), piece in board.items():
        if piece == 'O':  # White stone
            liberties = get_liberties(board, x, y)
            if len(liberties) == 1:  # Stone has only one liberty
                # Check if placing a black stone here would capture
                liberty_x, liberty_y = liberties[0]
                potential_captures.append((liberty_x, liberty_y))
    
    # Convert coordinates to Go notation
    go_moves = []
    cols = 'ABCDEFGHIJKL'
    for x, y in potential_captures:
        move = f"{cols[x]}{y+1}"
        go_moves.append(move)
    
    print("Potential capturing moves:", go_moves)

find_captures()