def create_board():
    board = [['.'] * 10 for _ in range(10)]
    # Add black stones (X)
    black_stones = [(9,3), (8,2), (8,4), (7,2), (7,4), (7,5), (6,3), (2,1), (2,4)]
    # Add white stones (O)
    white_stones = [(8,2), (8,3), (7,3), (7,6), (7,7), (4,2), (4,8), (2,4), (2,1)]
    
    for x, y in black_stones:
        board[9-x][y-1] = 'X'
    for x, y in white_stones:
        board[9-x][y-1] = 'O'
    return board

def get_liberties(board, x, y, visited):
    if x < 0 or x >= 10 or y < 0 or y >= 10:
        return set()
    if (x,y) in visited:
        return set()
    if board[x][y] == '.':
        return {(x,y)}
    
    visited.add((x,y))
    stone = board[x][y]
    liberties = set()
    
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 10 and 0 <= ny < 10:
            if board[nx][ny] == '.':
                liberties.add((nx,ny))
            elif board[nx][ny] == stone:
                liberties.update(get_liberties(board, nx, ny, visited))
    
    return liberties

def find_capturing_moves(board):
    potential_moves = []
    for x in range(10):
        for y in range(10):
            if board[x][y] == 'O':  # Look at white stones
                liberties = get_liberties(board, x, y, set())
                if len(liberties) == 1:  # If only one liberty
                    liberty = list(liberties)[0]
                    potential_moves.append((liberty, len(liberties)))
    
    return potential_moves

board = create_board()
moves = find_capturing_moves(board)

# Convert coordinates to Go notation
def coord_to_notation(x, y):
    return f"{chr(65+y)}{10-x}"

for move, liberties in moves:
    print(f"Potential capturing move at {coord_to_notation(move[0], move[1])}")