def is_valid_move(x, y, board):
    return 0 <= x < 4 and 0 <= y < 3 and board[y][x] == '.'

def knight_moves(x, y):
    return [(x+2, y+1), (x+2, y-1), (x-2, y+1), (x-2, y-1),
            (x+1, y+2), (x+1, y-2), (x-1, y+2), (x-1, y-2)]

def swap_knights(board, moves, turn):
    # Check if the current board matches the target configuration
    target_board = [['w', '.', '.', 'w'], ['.', '.', '.', '.'], ['B', 'B', '.', '.']]
    if board == target_board:
        return True, moves

    # Find current knight positions
    knight_positions = {'w': [], 'B': []}
    for y in range(3):
        for x in range(4):
            if board[y][x] in knight_positions:
                knight_positions[board[y][x]].append((x, y))

    current_knights = knight_positions[turn]

    for i, (x, y) in enumerate(current_knights):
        for nx, ny in knight_moves(x, y):
            if is_valid_move(nx, ny, board):
                # Make the move
                board[y][x], board[ny][nx] = '.', turn
                current_knights[i] = (nx, ny)
                moves.append(f"{turn},{chr(65+x)}{3-y},{chr(65+nx)}{3-ny}")

                # Alternate turn
                next_turn = 'B' if turn == 'w' else 'w'
                success, result_moves = swap_knights(board, moves, next_turn)
                if success:
                    return True, result_moves

                # Undo the move
                board[ny][nx], board[y][x] = '.', turn
                current_knights[i] = (x, y)
                moves.pop()

    return False, []

# Initial board setup
board = [['B', '.', '.', 'B'], ['.', '.', '.', '.'], ['w', 'w', '.', '.']]
success, result_moves = swap_knights(board, [], 'w')

if success:
    print(result_moves)
else:
    print("No")