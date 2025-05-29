def get_liberties(board, x, y, color, visited=None):
    if visited is None:
        visited = set()
    
    if x < 0 or x >= 12 or y < 0 or y >= 12:
        return set()
    
    if board[x][y] != color:
        if board[x][y] == '.':
            return {(x, y)}
        return set()
    
    if (x, y) in visited:
        return set()
    
    visited.add((x, y))
    liberties = set()
    
    # Check all adjacent positions
    for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
        new_x, new_y = x + dx, y + dy
        liberties.update(get_liberties(board, new_x, new_y, color, visited))
    
    return liberties

def create_board():
    board = [['.'] * 12 for _ in range(12)]
    # Input the board state (coordinates are reversed in array representation)
    stones = {
        'X': [(11,7), (10,6), (10,8), (8,2), (2,3), (1,2)],
        'O': [(11,2), (9,7), (9,4), (7,6), (4,0)]
    }
    # Add the upper right cluster
    white_cluster = [(11,7), (10,7), (10,6), (9,7)]
    black_cluster = [(11,6), (11,8), (10,5), (9,6), (9,8), (8,7)]
    
    for x, y in white_cluster:
        board[x][y] = 'O'
    for x, y in black_cluster:
        board[x][y] = 'X'
    for color in stones:
        for x, y in stones[color]:
            board[x][y] = color
            
    return board

def check_capture_moves():
    board = create_board()
    potential_moves = []
    
    # Check each empty point
    for i in range(12):
        for j in range(12):
            if board[i][j] == '.':
                # Try move here and check if it captures any white groups
                board[i][j] = 'X'
                # Check adjacent white groups
                for di, dj in [(0,1), (0,-1), (1,0), (-1,0)]:
                    ni, nj = i + di, j + dj
                    if 0 <= ni < 12 and 0 <= nj < 12 and board[ni][nj] == 'O':
                        liberties = get_liberties(board, ni, nj, 'O')
                        if not liberties:
                            # This move captures white stones
                            # Convert to board coordinates (A-L, 1-12)
                            move = f"{chr(ord('A') + j)}{12-i}"
                            potential_moves.append(move)
                board[i][j] = '.'
    
    print(f"Potential capturing moves: {potential_moves}")

check_capture_moves()