def is_valid_move(x, y, board):
    return 0 <= x < 4 and 0 <= y < 3 and board[y][x] == '.'

def get_knight_moves(x, y):
    moves = [(x+2, y+1), (x+2, y-1), (x-2, y+1), (x-2, y-1),
             (x+1, y+2), (x+1, y-2), (x-1, y+2), (x-1, y-2)]
    return moves

def swap_knights(board, moves, turn, max_depth=10):
    if board == [['.', 'w', '.', 'B'], ['.', '.', '.', 'B'], ['.', 'w', '.', '.']]:
        return True, moves

    if turn > max_depth:
        return False, []

    knight_positions = {'B': [], 'w': []}
    for y in range(3):
        for x in range(4):
            if board[y][x] in knight_positions:
                knight_positions[board[y][x]].append((x, y))

    current_color = 'B' if turn % 2 == 0 else 'w'
    for x, y in knight_positions[current_color]:
        for nx, ny in get_knight_moves(x, y):
            if is_valid_move(nx, ny, board):
                board[y][x], board[ny][nx] = '.', current_color
                moves.append(f"{current_color},{chr(65+x)}{3-y},{chr(65+nx)}{3-ny}")
                found, result = swap_knights(board, moves, turn + 1, max_depth)
                if found:
                    return True, result
                moves.pop()
                board[ny][nx], board[y][x] = '.', current_color

    return False, []

initial_board = [['.', 'B', '.', '.'], ['.', '.', '.', 'w'], ['.', 'B', '.', 'w']]
found, result = swap_knights(initial_board, [], 0)
if found:
    print(result)
else:
    print("No")