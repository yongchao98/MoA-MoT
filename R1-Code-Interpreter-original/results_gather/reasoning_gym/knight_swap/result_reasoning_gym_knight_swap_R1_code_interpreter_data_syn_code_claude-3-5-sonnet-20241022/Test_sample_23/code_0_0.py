def get_valid_moves(pos, board):
    moves = []
    x, y = pos
    knight_moves = [(-2,-1),(-2,1),(-1,-2),(-1,2),(1,-2),(1,2),(2,-1),(2,1)]
    for dx, dy in knight_moves:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < 4 and 0 <= new_y < 4 and board[new_y][new_x] == '.':
            moves.append((new_x, new_y))
    return moves

def is_complete(board, target_w, target_b):
    w_positions = [(x,y) for y in range(4) for x in range(4) if board[y][x] == 'w']
    b_positions = [(x,y) for y in range(4) for x in range(4) if board[y][x] == 'B']
    return sorted(w_positions) == sorted(target_b) and sorted(b_positions) == sorted(target_w)

def solve(board, moves, white_turn, target_w, target_b, depth=0):
    if depth > 20:  # Prevent infinite recursion
        return None
    if is_complete(board, target_w, target_b):
        return moves
    
    pieces = [(x,y) for y in range(4) for x in range(4) if board[y][x] == ('w' if white_turn else 'B')]
    for piece in pieces:
        valid_moves = get_valid_moves(piece, board)
        for move in valid_moves:
            # Make move
            new_board = [row[:] for row in board]
            new_board[piece[1]][piece[0]] = '.'
            new_board[move[1]][move[0]] = 'w' if white_turn else 'B'
            
            # Record move
            move_str = f"{'w' if white_turn else 'B'},{chr(65+piece[0])}{piece[1]+1},{chr(65+move[0])}{move[1]+1}"
            
            # Recursive solve
            result = solve(new_board, moves + [move_str], not white_turn, target_w, target_b, depth+1)
            if result:
                return result
    return None

# Initial board
board = [
    ['w', 'w', '.', '.'],
    ['.', '.', '.', 'B'],
    ['.', 'B', '.', '.']
]

# Target positions
white_start = [(0,3), (1,3)]  # A3, B3
black_start = [(1,1), (3,2)]  # B1, D2

# Find solution
solution = solve(board, [], True, white_start, black_start)
if solution:
    print(solution)
else:
    print("No")