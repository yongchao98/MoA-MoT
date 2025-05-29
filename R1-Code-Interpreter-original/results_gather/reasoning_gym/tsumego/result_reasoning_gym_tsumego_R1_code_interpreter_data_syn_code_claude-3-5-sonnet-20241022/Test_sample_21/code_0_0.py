def create_board():
    board = {}
    # Initialize the given position
    positions = {
        'O': ['C13', 'H10', 'H8', 'H2', 'K11', 'L11', 'M11', 'L10', 'L6', 'M5'],
        'X': ['F13', 'L13', 'C10', 'D8', 'K10', 'M10', 'L9', 'L8', 'X2']
    }
    
    for color in positions:
        for pos in positions[color]:
            col = ord(pos[0]) - ord('A')
            row = int(pos[1:]) - 1
            board[(col, row)] = color
    
    return board

def get_liberties(board, x, y):
    liberties = []
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < 13 and 0 <= new_y < 13:
            if (new_x, new_y) not in board:
                liberties.append((new_x, new_y))
    return liberties

def analyze_position():
    board = create_board()
    # Analyze white group in upper right
    white_group = [(10,10), (11,10), (12,10)]  # K11, L11, M11
    all_liberties = set()
    
    for stone in white_group:
        libs = get_liberties(board, stone[0], stone[1])
        for lib in libs:
            all_liberties.add(lib)
    
    # Convert coordinates back to board notation
    liberty_points = []
    for lib in all_liberties:
        col = chr(lib[0] + ord('A'))
        row = str(lib[1] + 1)
        liberty_points.append(col + row)
    
    print("White group liberties:", sorted(liberty_points))

analyze_position()