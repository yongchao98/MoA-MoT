def create_board():
    board = {}
    board_str = """
   A B C D E F G H I J K
11 . . X . . . . X . X .
10 X . . . X . . . X O X
 9 . . . . . . . X O O .
 8 . . . . . . . . X O X
 7 . . . . . . . . X X .
 6 . . . . . . . . . X .
 5 . . X . . . . . . . .
 4 . . . . . . . O . . .
 3 . . . . . . . . . . X
 2 . . . . . . . . . . .
 1 . . . . . O . X . . .
"""
    lines = board_str.strip().split('\n')
    for i, line in enumerate(lines[1:], 1):
        row = 12 - i
        cols = line[1:].split()
        for j, piece in enumerate(cols):
            if piece != '.':
                col = chr(ord('A') + j)
                board[(col, row)] = piece
    return board

def get_liberties(board, pos):
    col, row = pos
    liberties = []
    for dc, dr in [(0,1), (0,-1), (1,0), (-1,0)]:
        new_col = chr(ord(col) + dc)
        new_row = row + dr
        if 'A' <= new_col <= 'K' and 1 <= new_row <= 11:
            if (new_col, new_row) not in board:
                liberties.append((new_col, new_row))
    return liberties

def analyze_white_groups(board):
    white_stones = [(pos, stone) for pos, stone in board.items() if stone == 'O']
    for pos, _ in white_stones:
        liberties = get_liberties(board, pos)
        if len(liberties) == 1:
            print(f"White stone at {pos} has only one liberty at {liberties[0]}")

board = create_board()
analyze_white_groups(board)