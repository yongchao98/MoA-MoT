def create_board():
    board = [['.'] * 11 for _ in range(11)]
    # Fill known positions
    black_stones = [(0,0), (1,0), (6,8), (7,8), (5,8), (5,6), (7,6), (5,3), (5,1), (8,4), (10,6)]
    white_stones = [(1,7), (3,6), (6,7), (6,8), (6,6), (8,7), (8,6), (0,3)]
    
    for x, y in black_stones:
        board[10-y][x] = 'X'
    for x, y in white_stones:
        board[10-y][x] = 'O'
    return board

def get_liberties(board, x, y):
    liberties = []
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < 11 and 0 <= new_y < 11 and board[new_y][new_x] == '.':
            liberties.append((new_x, new_y))
    return liberties

def analyze_position():
    board = create_board()
    white_group = [(6,7), (6,8), (6,6)]  # Core white stones
    all_liberties = set()
    
    for x, y in white_group:
        libs = get_liberties(board, x, 10-y)
        for lib in libs:
            all_liberties.add(lib)
    
    print("Critical points to consider:")
    for x, y in all_liberties:
        print(f"({chr(65+x)},{11-y})")

analyze_position()