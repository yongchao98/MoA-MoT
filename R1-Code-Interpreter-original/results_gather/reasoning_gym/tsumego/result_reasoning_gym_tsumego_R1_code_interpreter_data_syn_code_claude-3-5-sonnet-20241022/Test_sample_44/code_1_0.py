def create_board():
    # Create empty 11x11 board
    board = [['.'] * 11 for _ in range(11)]
    
    # Place Black stones (X)
    black_positions = [
        (4,10), (6,10),  # row 10
        (3,9), (7,9),    # row 9
        (4,8), (6,8),    # row 8
        (0,7),           # row 7
        (4,6),           # row 6
        (3,5),           # row 5
        (0,4), (6,4), (7,4),  # row 4
        (0,3),           # row 3
        (6,1)            # row 1
    ]
    
    # Place White stones (O)
    white_positions = [
        (5,11),          # row 11
        (4,9), (5,9), (6,9),  # row 9
        (5,8), (10,8),   # row 8
        (9,5)            # row 5
    ]
    
    # Place stones on board
    for x, y in black_positions:
        board[11-y][x] = 'X'
    for x, y in white_positions:
        board[11-y][x] = 'O'
    
    return board

def print_board(board):
    print("   A B C D E F G H I J K")
    for i in range(11):
        row_num = 11 - i
        print(f"{row_num:2d} {' '.join(board[i])}")

def check_capture_moves(board):
    captures = []
    directions = [(0,1), (0,-1), (1,0), (-1,0)]
    
    def count_group_liberties(x, y, visited=None):
        if visited is None:
            visited = set()
        
        if x < 0 or x >= 11 or y < 0 or y >= 11:
            return set()
            
        if board[y][x] == '.':
            return {(x,y)}
            
        if (x,y) in visited or board[y][x] != 'O':
            return set()
            
        visited.add((x,y))
        liberties = set()
        
        for dx, dy in directions:
            new_x, new_y = x + dx, y + dy
            liberties.update(count_group_liberties(new_x, new_y, visited))
        
        return liberties

    # Check each position
    for y in range(11):
        for x in range(11):
            if board[y][x] == 'O':
                liberties = count_group_liberties(x, y)
                if len(liberties) <= 2:  # Vulnerable group
                    for liberty in liberties:
                        # Check if this liberty point could be a capturing move
                        captures.append((liberty[0], liberty[1], len(liberties)))

    return captures

board = create_board()
print_board(board)
print("\nPotential capturing moves:")
captures = check_capture_moves(board)
for x, y, liberties in captures:
    col = chr(ord('A') + x)
    row = 11 - y
    print(f"Move at {col}{row} (liberties: {liberties})")