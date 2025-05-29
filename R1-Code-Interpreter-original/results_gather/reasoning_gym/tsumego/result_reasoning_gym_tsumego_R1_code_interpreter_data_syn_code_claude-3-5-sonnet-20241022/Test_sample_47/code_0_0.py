def create_board():
    board = [['.'] * 10 for _ in range(10)]
    # Place black stones (X)
    black_positions = [(9,3), (8,3), (8,6), (9,7), (5,5), (5,8), (5,9), 
                      (4,7), (3,2), (3,4), (3,7), (2,8)]
    # Place white stones (O)
    white_positions = [(6,3), (4,7), (3,0), (3,5), (3,8), (4,8), (4,9)]
    
    for x, y in black_positions:
        board[x][y] = 'X'
    for x, y in white_positions:
        board[x][y] = 'O'
    return board

def get_liberties(board, x, y):
    liberties = 0
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < 10 and 0 <= new_y < 10 and board[new_x][new_y] == '.':
            liberties += 1
    return liberties

def find_vulnerable_groups(board):
    vulnerable = []
    for i in range(10):
        for j in range(10):
            if board[i][j] == 'O':  # Check white stones
                libs = get_liberties(board, i, j)
                if libs <= 2:  # Stones with 2 or fewer liberties are vulnerable
                    vulnerable.append((i, j, libs))
    return vulnerable

def find_capturing_moves(board):
    potential_moves = []
    vulnerable = find_vulnerable_groups(board)
    
    # Check all empty points around vulnerable stones
    for x, y, libs in vulnerable:
        for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
            new_x, new_y = x + dx, y + dy
            if 0 <= new_x < 10 and 0 <= new_y < 10 and board[new_x][new_y] == '.':
                # Convert to board coordinates (A-J, 1-10)
                move = (chr(ord('A') + new_y), 10 - new_x)
                potential_moves.append((move, libs))
    
    return potential_moves

board = create_board()
moves = find_capturing_moves(board)
print("Potential capturing moves:", moves)