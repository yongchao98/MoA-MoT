def create_board():
    board = [['.'] * 9 for _ in range(9)]
    # Place black stones (X)
    black_stones = [(8,4), (8,0), (7,2), (5,3), (5,4), (5,7), 
                    (4,1), (3,2), (3,4), (2,1), (2,3)]
    # Place white stones (O)
    white_stones = [(5,6), (4,2), (4,3), (4,5), (3,3), (3,7),
                    (1,3), (1,5)]
    
    for x, y in black_stones:
        board[x][y] = 'X'
    for x, y in white_stones:
        board[x][y] = 'O'
    return board

def get_liberties(board, x, y, visited=None):
    if visited is None:
        visited = set()
    
    if not (0 <= x < 9 and 0 <= y < 9):
        return set()
    
    if board[x][y] == '.':
        return {(x,y)}
    
    if (x,y) in visited:
        return set()
    
    visited.add((x,y))
    stone = board[x][y]
    liberties = set()
    
    for dx, dy in [(0,1), (1,0), (0,-1), (-1,0)]:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 9 and 0 <= ny < 9:
            if board[nx][ny] == '.':
                liberties.add((nx,ny))
            elif board[nx][ny] == stone:
                liberties.update(get_liberties(board, nx, ny, visited))
    
    return liberties

def find_capturing_moves():
    board = create_board()
    potential_moves = []
    
    # Check each white stone
    for x in range(9):
        for y in range(9):
            if board[x][y] == 'O':
                liberties = get_liberties(board, x, y)
                if len(liberties) == 1:
                    liberty = liberties.pop()
                    potential_moves.append((liberty, (x,y)))
    
    # Sort moves by potential capture size
    moves_with_score = []
    for move, target in potential_moves:
        # Simple scoring: count connected stones that might be captured
        score = 1  # At least one stone
        moves_with_score.append((move, score))
    
    moves_with_score.sort(key=lambda x: x[1], reverse=True)
    
    # Convert to board coordinates
    if moves_with_score:
        x, y = moves_with_score[0][0]
        col = chr(ord('A') + y)
        row = 9 - x
        print(f"{col}{row}")

find_capturing_moves()