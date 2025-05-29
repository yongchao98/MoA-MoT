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

def check_capture_size(board, x, y, visited=None):
    if visited is None:
        visited = set()
    
    if not (0 <= x < 10 and 0 <= y < 10) or board[x][y] != 'O':
        return 0
    
    if (x, y) in visited:
        return 0
        
    visited.add((x, y))
    count = 1
    
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        new_x, new_y = x + dx, y + dy
        count += check_capture_size(board, new_x, new_y, visited)
    
    return count

def analyze_moves(board):
    moves = {}
    for i in range(10):
        for j in range(10):
            if board[i][j] == '.':
                capture_potential = 0
                for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
                    new_x, new_y = i + dx, j + dy
                    if 0 <= new_x < 10 and 0 <= new_y < 10 and board[new_x][new_y] == 'O':
                        liberties = get_liberties(board, new_x, new_y)
                        if liberties == 1:
                            capture_potential += check_capture_size(board, new_x, new_y)
                
                if capture_potential > 0:
                    move = (chr(ord('A') + j), 10 - i)
                    moves[move] = capture_potential

    return moves

# Create board and analyze
board = create_board()
capturing_moves = analyze_moves(board)
print("Moves and number of stones they can capture:", capturing_moves)

# Additional analysis: Check if moves are adjacent to existing black stones
def check_black_connection(board, x, y):
    count = 0
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < 10 and 0 <= new_y < 10 and board[new_x][new_y] == 'X':
            count += 1
    return count

# Check connections for each capturing move
for i in range(10):
    for j in range(10):
        if board[i][j] == '.':
            connections = check_black_connection(board, i, j)
            if connections > 0:
                move = (chr(ord('A') + j), 10 - i)
                print(f"Move {move} connects to {connections} black stones")