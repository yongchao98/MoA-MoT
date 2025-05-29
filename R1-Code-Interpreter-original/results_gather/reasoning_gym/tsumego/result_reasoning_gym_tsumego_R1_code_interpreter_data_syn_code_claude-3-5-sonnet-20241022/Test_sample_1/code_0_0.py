def create_board():
    board = {}
    board_str = """
   A B C D E F G H I J K L M
13 . . . . . . . . . . . O .
12 . . . . . . . . . . . . .
11 . . . . . . . . . . X . .
10 . . . . . . . . . . . . .
 9 . . . . . O . O X . . . .
 8 . . . . . . . . . . O . .
 7 . . . . . . . . . . . . .
 6 . . . . O . . . . . . . .
 5 . . . . . . . . X . . . .
 4 . . . . X . . . . . . . .
 3 . X . X O X . . . O . . .
 2 . . X O O . X . . . X . .
 1 . . . X O X . . . . . . .
"""
    rows = board_str.strip().split('\n')
    for i, row in enumerate(rows[1:], 1):  # Skip header row
        cols = row.split()[1:]  # Skip row number
        for j, stone in enumerate(cols):
            if stone != '.':
                col = chr(ord('A') + j)
                row_num = 14 - i  # Convert to actual row numbers
                board[f"{col}{row_num}"] = stone
    return board

def get_adjacent(pos):
    col, row = pos[0], int(pos[1:])
    adj = []
    # Check all four directions
    if col > 'A': adj.append(f"{chr(ord(col)-1)}{row}")
    if col < 'M': adj.append(f"{chr(ord(col)+1)}{row}")
    if row > 1: adj.append(f"{col}{row-1}")
    if row < 13: adj.append(f"{col}{row+1}")
    return adj

def find_capturing_moves():
    board = create_board()
    potential_moves = []
    
    # Check each white stone's liberties
    for pos, stone in board.items():
        if stone == 'O':
            liberties = 0
            empty_adj = []
            for adj in get_adjacent(pos):
                if adj not in board:
                    liberties += 1
                    empty_adj.append(adj)
            if liberties == 1:  # If stone has only one liberty
                potential_moves.extend(empty_adj)
    
    return potential_moves

# Find and print potential capturing moves
moves = find_capturing_moves()
print("Potential capturing moves:", moves)