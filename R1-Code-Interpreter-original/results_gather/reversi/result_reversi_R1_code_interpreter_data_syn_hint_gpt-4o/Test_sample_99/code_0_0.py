def print_board(board):
    for row in board:
        print(''.join(row))

def is_on_board(x, y, n):
    return 0 <= x < n and 0 <= y < n

def flip_pieces(board, x, y, player, opponent, n):
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
    to_flip = []
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        pieces_to_flip = []
        while is_on_board(nx, ny, n) and board[nx][ny] == opponent:
            pieces_to_flip.append((nx, ny))
            nx += dx
            ny += dy
        if is_on_board(nx, ny, n) and board[nx][ny] == player:
            to_flip.extend(pieces_to_flip)
    for fx, fy in to_flip:
        board[fx][fy] = player

def play_game():
    n = 6
    board = [
        ['*', '*', '*', '*', '*', '*'],
        ['*', '1', '0', '*', '*', '*'],
        ['*', '0', '0', '0', '*', '*'],
        ['*', '*', '0', '1', '*', '*'],
        ['*', '*', '*', '*', '*', '*'],
        ['*', '*', '*', '*', '*', '*']
    ]

    # Round 1: Your's landing position is (2,3)
    board[2][3] = '0'
    flip_pieces(board, 2, 3, '0', '1', n)

    # Round 2: The opponent's landing position is (2,2)
    board[2][2] = '1'
    flip_pieces(board, 2, 2, '1', '0', n)

    # Round 3: Your's landing position is (3,2)
    board[3][2] = '0'
    flip_pieces(board, 3, 2, '0', '1', n)

    # Print the final board
    print_board(board)

play_game()