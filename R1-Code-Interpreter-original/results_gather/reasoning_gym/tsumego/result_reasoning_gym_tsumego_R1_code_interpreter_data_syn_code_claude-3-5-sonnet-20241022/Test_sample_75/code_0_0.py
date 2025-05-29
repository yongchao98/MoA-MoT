def create_board():
    board = [['.'] * 10 for _ in range(10)]
    # Add black stones (X)
    black_stones = [(8,2), (8,3), (8,4), (5,0), (5,1), (5,2), (4,0), (4,3), (3,0), (2,1)]
    # Add white stones (O)
    white_stones = [(9,4), (7,2), (5,1), (5,2), (4,1), (3,7), (5,6), (3,7), (0,9)]
    
    for x, y in black_stones:
        board[x][y] = 'X'
    for x, y in white_stones:
        board[x][y] = 'O'
    
    return board

def check_liberties(board, x, y):
    if x < 0 or x >= 10 or y < 0 or y >= 10:
        return []
    liberties = []
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < 10 and 0 <= new_y < 10 and board[new_x][new_y] == '.':
            liberties.append((new_x, new_y))
    return liberties

def analyze_position():
    board = create_board()
    # Check liberties of white stones at B5 and C5
    white_group_liberties = set()
    for x, y in [(5,1), (5,2)]:  # B5 and C5 coordinates
        white_group_liberties.update(check_liberties(board, x, y))
    
    # Convert coordinates to Go notation
    for x, y in white_group_liberties:
        col = chr(y + ord('A'))
        row = 10 - x
        print(f"Liberty found at: {col}{row}")

analyze_position()