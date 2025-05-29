def create_board():
    board = [['.'] * 10 for _ in range(10)]
    # Add black stones (X)
    black_stones = [(8,2), (7,3), (6,3), (5,4), (4,1), (3,1), (2,2), (2,8), (8,0)]
    # Add white stones (O)
    white_stones = [(7,1), (7,2), (5,1), (5,2), (5,0), (6,0), (8,5)]
    
    for x, y in black_stones:
        board[y][x] = 'X'
    for x, y in white_stones:
        board[y][x] = 'O'
    return board

def get_liberties(board, x, y, visited=None):
    if visited is None:
        visited = set()
    
    if x < 0 or x >= 10 or y < 0 or y >= 10:
        return set()
    
    stone = board[y][x]
    if stone == '.':
        return {(x,y)}
    
    if stone not in 'OX' or (x,y) in visited:
        return set()
    
    visited.add((x,y))
    liberties = set()
    
    # Check all adjacent positions
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        new_x, new_y = x + dx, y + dy
        liberties.update(get_liberties(board, new_x, new_y, visited))
    
    return liberties

def find_capturing_moves(board):
    potential_moves = []
    
    # Find all white groups and their liberties
    visited = set()
    for y in range(10):
        for x in range(10):
            if board[y][x] == 'O' and (x,y) not in visited:
                liberties = get_liberties(board, x, y)
                if len(liberties) == 1:
                    liberty = liberties.pop()
                    potential_moves.append((liberty, (x,y)))
    
    return potential_moves

board = create_board()
moves = find_capturing_moves(board)

# Convert coordinates to Go notation
def coord_to_notation(x, y):
    return f"{chr(65+x)}{10-y}"

for move, stone in moves:
    print(f"Move at {coord_to_notation(move[0], move[1])} can capture stone at {coord_to_notation(stone[0], stone[1])}")