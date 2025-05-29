def print_board(board):
    print("   A B C D E F G H I J K")
    for i in range(11, 0, -1):
        row = f"{i:2d} "
        for j in range(11):
            row += board[i-1][j] + " "
        print(row)

def is_valid_move(x, y):
    return 0 <= x < 11 and 0 <= y < 11

def find_vital_points():
    board = [
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', 'X', '.', '.', 'X', '.'],
        ['.', '.', '.', 'O', '.', '.', '.', 'X', 'X', 'O', 'X'],
        ['.', '.', 'O', '.', '.', '.', '.', 'X', 'O', 'O', '.'],
        ['.', '.', '.', '.', '.', '.', 'X', '.', 'X', 'O', 'X'],
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.'],
        ['.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', 'O', '.', '.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.'],
        ['.', '.', '.', '.', '.', '.', '.', '.', 'X', 'O', '.']
    ]

    vital_points = []
    
    # Look for cutting points and vital points
    for y in range(11):
        for x in range(11):
            if board[y][x] == '.':
                # Count adjacent stones
                black_adj = 0
                white_adj = 0
                diag_black = 0
                diag_white = 0
                
                # Check orthogonal adjacencies
                for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
                    nx, ny = x+dx, y+dy
                    if is_valid_move(nx, ny):
                        if board[ny][nx] == 'X':
                            black_adj += 1
                        elif board[ny][nx] == 'O':
                            white_adj += 1
                
                # Check diagonal adjacencies
                for dx, dy in [(1,1), (1,-1), (-1,1), (-1,-1)]:
                    nx, ny = x+dx, y+dy
                    if is_valid_move(nx, ny):
                        if board[ny][nx] == 'X':
                            diag_black += 1
                        elif board[ny][nx] == 'O':
                            diag_white += 1
                
                # Identify potential vital points
                if (black_adj >= 1 and white_adj >= 1) or \
                   (black_adj >= 2) or \
                   (white_adj >= 2) or \
                   (black_adj + white_adj >= 2 and diag_black + diag_white >= 1):
                    col = chr(ord('A') + x)
                    row = 11 - y
                    score = black_adj * 2 + white_adj * 2 + diag_black + diag_white
                    vital_points.append((f"{col}{row}", score, 
                                      f"B:{black_adj} W:{white_adj} DB:{diag_black} DW:{diag_white}"))

    # Sort vital points by their strategic score
    vital_points.sort(key=lambda x: x[1], reverse=True)
    for point in vital_points:
        print(f"Point: {point[0]}, Score: {point[1]}, Adjacent stones: {point[2]}")

find_vital_points()