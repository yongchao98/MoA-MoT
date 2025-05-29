def create_board():
    board = [['.'] * 13 for _ in range(13)]
    # Add the stones (coordinates are 0-based)
    # White stones (O)
    white_positions = [(1,6), (1,5), (6,6), (1,7), (8,6), (9,12)]
    # Black stones (X)
    black_positions = [(0,7), (1,8), (2,7), (2,5), (0,5), (2,4), (3,6), 
                      (6,6), (10,1), (10,2), (10,10), (10,11), (4,11), (12,7)]
    
    for x, y in white_positions:
        board[y][x] = 'O'
    for x, y in black_positions:
        board[y][x] = 'X'
    
    return board

def get_liberties(board, x, y):
    liberties = []
    directions = [(0,1), (0,-1), (1,0), (-1,0)]
    for dx, dy in directions:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < 13 and 0 <= new_y < 13:
            if board[new_y][new_x] == '.':
                liberties.append((new_x, new_y))
    return liberties

def analyze_position():
    board = create_board()
    
    # Check each white stone for vulnerabilities
    for y in range(13):
        for x in range(13):
            if board[y][x] == 'O':
                liberties = get_liberties(board, x, y)
                if len(liberties) == 1:  # If stone has only one liberty
                    liberty = liberties[0]
                    print(f"White stone at ({x+1},{13-y}) has only one liberty at {chr(65+liberty[0])}{13-liberty[1]}")

analyze_position()