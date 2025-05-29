def create_board():
    board = [['.'] * 10 for _ in range(10)]
    # Place black stones (X)
    black_stones = [(9,2), (8,0), (7,1), (7,3), (7,9), (6,2), (6,9), (5,9), (4,9), (1,1), (1,8)]
    # Place white stones (O)
    white_stones = [(8,1), (8,2), (7,2), (5,5), (5,8), (3,3), (2,8)]
    
    for x, y in black_stones:
        board[x][y] = 'X'
    for x, y in white_stones:
        board[x][y] = 'O'
    return board

def get_liberties(board, x, y, visited=None):
    if visited is None:
        visited = set()
    
    if x < 0 or x >= 10 or y < 0 or y >= 10:
        return set()
    
    if board[x][y] == '.':
        return {(x,y)}
    
    if (x,y) in visited:
        return set()
    
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

def find_capturing_moves():
    board = create_board()
    best_moves = []
    max_captures = 0
    
    # Find white groups and their liberties
    for i in range(10):
        for j in range(10):
            if board[i][j] == 'O':
                liberties = get_liberties(board, i, j)
                if len(liberties) <= 2:  # Vulnerable group
                    for liberty in liberties:
                        # Count stones that would be captured
                        captured = 1  # At least one stone
                        for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
                            nx, ny = i + dx, j + dy
                            if 0 <= nx < 10 and 0 <= ny < 10 and board[nx][ny] == 'O':
                                captured += 1
                        if captured >= max_captures:
                            max_captures = captured
                            best_moves.append((liberty, captured))
    
    # Convert to board coordinates
    result = []
    for (x,y), captures in best_moves:
        if captures == max_captures:
            col = chr(y + ord('A'))
            row = 10 - x
            result.append(f"{col}{row}")
    
    print(f"Best moves: {result} capturing {max_captures} stones")

find_capturing_moves()