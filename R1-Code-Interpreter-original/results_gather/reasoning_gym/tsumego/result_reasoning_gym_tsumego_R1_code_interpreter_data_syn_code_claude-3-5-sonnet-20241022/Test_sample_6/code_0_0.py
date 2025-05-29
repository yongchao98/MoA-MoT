def create_board():
    board = {}
    board_str = """   A B C D E F G H I J K L
12 . . . . . . X . . O . .
11 . . . . . . . . . . . .
10 . . . . . . O . . . . .
 9 . . . . . . . . . . . .
 8 . . . . . . . O . . . .
 7 . . . . . . . . . . . X
 6 . . . . . . . . . . . X
 5 . . . . . . . . . . . .
 4 . . . . . . . . . X O .
 3 . O . . . X X . X O X .
 2 . . . . . . . X O O . X
 1 . . . . . O . . X O X ."""
    
    rows = board_str.strip().split('\n')[1:]  # Skip the column labels
    for row_num, row in enumerate(rows):
        row_pieces = row[3:].split()  # Skip the row number
        for col_num, piece in enumerate(row_pieces):
            if piece != '.':
                board[(col_num, 11-row_num)] = piece
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
        if piece == 'O':
            liberties = get_liberties(board, x, y)
            if len(liberties) == 1:  # Stone has only one liberty
                liberty_x, liberty_y = liberties[0]
                # Check if placing a black stone here would capture
                potential_captures.append((liberty_x, liberty_y))
    
    # Convert coordinates to Go notation
    for x, y in potential_captures:
        col = chr(65 + x)
        row = y + 1
        print(f"Potential capturing move: {col}{row}")

find_captures()