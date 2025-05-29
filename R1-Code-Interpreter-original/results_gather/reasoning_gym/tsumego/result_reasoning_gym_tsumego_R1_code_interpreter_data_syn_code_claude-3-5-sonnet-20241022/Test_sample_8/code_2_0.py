def check_capture(board_str):
    # Convert the board string into a 2D list for easier processing
    lines = [line.strip() for line in board_str.split('\n') if line.strip() and not line.strip().startswith('   ')]
    board = []
    for line in lines:
        row = line[4:].split()  # Skip the row number and split
        board.append(row)
    
    # Function to check if a position is on the board
    def is_valid(x, y):
        return 0 <= x < 12 and 0 <= y < 12
    
    # Function to get liberties of a stone
    def get_liberties(x, y):
        liberties = []
        for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
            new_x, new_y = x + dy, y + dx  # Note: swapped dx,dy because board is in row,col format
            if is_valid(new_x, new_y) and board[new_x][new_y] == '.':
                liberties.append((new_x, new_y))
        return liberties
    
    # Check each white stone for potential capture
    captures = []
    for i in range(12):
        for j in range(12):
            if board[i][j] == 'O':  # Found a white stone
                # Check if it's nearly surrounded by black stones
                surrounded = 0
                liberties = get_liberties(i, j)
                for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
                    new_x, new_y = i + dy, j + dx
                    if is_valid(new_x, new_y) and board[new_x][new_y] == 'X':
                        surrounded += 1
                
                if surrounded >= 3 and len(liberties) == 1:
                    # This stone is capturable
                    liberty = liberties[0]
                    # Convert to board coordinates (A-L, 12-1)
                    move = f"{chr(ord('A') + liberty[1])}{12-liberty[0]}"
                    captures.append((move, i, j))
    
    print("Potential capturing moves:", captures)

# The board state
board_str = """
   A B C D E F G H I J K L
12 . . O . . . . X . . . .
11 . . . . . . X O X . . .
10 . . . . . X O O . . . .
 9 . . . . . . X O X O . .
 8 . . . . . . . X . . O .
 7 . X O . O . . . . . . .
 6 . . . . . . . . . . . .
 5 O . . . . . . . . . . .
 4 . . . . . . . . . . . .
 3 . . . . . . . . X . . .
 2 . . . X . . . . O . . .
 1 . . X . . . . . . . . .
"""

check_capture(board_str)