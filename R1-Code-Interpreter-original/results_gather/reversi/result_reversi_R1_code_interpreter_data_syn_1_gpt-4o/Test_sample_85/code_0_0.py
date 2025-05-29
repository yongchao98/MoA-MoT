def print_board(board):
    for row in board:
        print(''.join(row))

def is_on_board(x, y, n):
    return 0 <= x < n and 0 <= y < n

def flip_pieces(board, start_x, start_y, color, n):
    opponent_color = '1' if color == '0' else '0'
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
    to_flip = []

    for dx, dy in directions:
        x, y = start_x, start_y
        x += dx
        y += dy
        pieces_to_flip = []

        while is_on_board(x, y, n) and board[x][y] == opponent_color:
            pieces_to_flip.append((x, y))
            x += dx
            y += dy

        if is_on_board(x, y, n) and board[x][y] == color:
            to_flip.extend(pieces_to_flip)

    for x, y in to_flip:
        board[x][y] = color

def play_game(n, moves):
    board = [
        ['*', '*', '*', '*', '*', '*'],
        ['*', '*', '0', '0', '0', '*'],
        ['*', '*', '0', '0', '*', '*'],
        ['*', '*', '0', '1', '*', '*'],
        ['*', '*', '*', '*', '*', '*'],
        ['*', '*', '*', '*', '*', '*']
    ]

    for move in moves:
        x, y, color = move
        board[x][y] = color
        flip_pieces(board, x, y, color, n)

    return board

n = 6
moves = [(1, 2, '0'), (1, 3, '1'), (1, 4, '0')]
final_board = play_game(n, moves)

# Flatten the board to the required output format
flattened_board = [cell for row in final_board for cell in row]
print(flattened_board)