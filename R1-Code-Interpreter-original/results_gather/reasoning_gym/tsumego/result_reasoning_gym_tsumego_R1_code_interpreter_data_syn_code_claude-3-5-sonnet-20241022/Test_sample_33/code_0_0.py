def create_board():
    board = {}
    board_str = """
   A B C D E F G H I J K L
12 . . . . . . . . . . . .
11 . . . . . . . . X . . .
10 O X X . . . . . O . O X
 9 X O O X . . . X . . X .
 8 X O . . O . . . . . . .
 7 . X . . . . . . . . . .
 6 . . . . . . . . . . . .
 5 . . . . . . . . O . . .
 4 . . . . . . . . . . . .
 3 . . . . . . . . . . X .
 2 . . . . . . O . . . . .
 1 X . . . . . . . . . . .
"""
    rows = board_str.strip().split('\n')[1:]  # Skip the first line
    for row_num, row in enumerate(rows):
        row_pieces = row[4:].split()  # Skip the row number
        for col_num, piece in enumerate(row_pieces):
            if piece != '.':
                board[(12 - row_num, chr(ord('A') + col_num))] = piece
    return board

def get_liberties(board, pos):
    row, col = pos
    liberties = []
    for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
        new_row = row + dr
        new_col = chr(ord(col) + dc)
        if 1 <= new_row <= 12 and 'A' <= new_col <= 'L':
            if (new_row, new_col) not in board:
                liberties.append((new_row, new_col))
    return liberties

def find_capturing_moves(board):
    captures = []
    for pos, piece in board.items():
        if piece == 'O':  # Look for white stones
            liberties = get_liberties(board, pos)
            if len(liberties) == 1:  # If only one liberty left
                # Check if this liberty is not shared with connected white stones
                # This is a simplified check
                liberty = liberties[0]
                captures.append(liberty)
    
    return captures

board = create_board()
potential_captures = find_capturing_moves(board)

for move in potential_captures:
    row, col = move
    print(f"Potential capturing move at {col}{row}")