def is_valid_move(x, y, board):
    return 0 <= x < 4 and 0 <= y < 3 and board[y][x] == '.'

def knight_moves(x, y):
    return [(x + 2, y + 1), (x + 2, y - 1), (x - 2, y + 1), (x - 2, y - 1),
            (x + 1, y + 2), (x + 1, y - 2), (x - 1, y + 2), (x - 1, y - 2)]

def swap_knights(board, moves, turn, depth=0):
    # Base case: check if the board is in the target configuration
    if board == [['.', 'B', 'B', '.'], ['.', '.', '.', '.'], ['w', '.', 'w', '.']]:
        return True, moves

    # Limit the depth to prevent infinite recursion
    if depth > 20:
        return False, []

    for y in range(3):
        for x in range(4):
            if (turn == 'w' and board[y][x] == 'w') or (turn == 'B' and board[y][x] == 'B'):
                for nx, ny in knight_moves(x, y):
                    if is_valid_move(nx, ny, board):
                        # Make the move
                        board[y][x], board[ny][nx] = '.', board[y][x]
                        moves.append(f"{turn},{chr(65+x)}{3-y},{chr(65+nx)}{3-ny}")
                        # Alternate turn
                        next_turn = 'B' if turn == 'w' else 'w'
                        solved, result_moves = swap_knights(board, moves, next_turn, depth + 1)
                        if solved:
                            return True, result_moves
                        # Backtrack
                        moves.pop()
                        board[ny][nx], board[y][x] = '.', board[ny][nx]
    return False, []

# Initial board setup
board = [['w', '.', 'w', '.'], ['.', '.', '.', '.'], ['.', 'B', 'B', '.']]
moves = []
solved, result_moves = swap_knights(board, moves, 'w')

if solved:
    print(result_moves)
else:
    print("No")