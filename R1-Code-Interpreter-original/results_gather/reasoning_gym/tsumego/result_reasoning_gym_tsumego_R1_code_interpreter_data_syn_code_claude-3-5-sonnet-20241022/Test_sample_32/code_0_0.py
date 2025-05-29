def create_board():
    board = [['.'] * 12 for _ in range(12)]
    # Place black stones (X)
    black_stones = [(11,1), (9,8), (8,8), (9,11), (9,9), (7,7), (5,3), (3,2), (6,7)]
    # Place white stones (O)
    white_stones = [(10,0), (10,10), (9,9), (8,9), (7,9), (2,0), (1,3), (10,8), (7,3)]
    
    for x, y in black_stones:
        board[x][y] = 'X'
    for x, y in white_stones:
        board[x][y] = 'O'
    return board

def count_liberties(board, x, y, visited=None):
    if visited is None:
        visited = set()
    
    if x < 0 or x >= 12 or y < 0 or y >= 12:
        return 0
    
    if board[x][y] == '.':
        return 1
    
    if (x, y) in visited:
        return 0
        
    visited.add((x, y))
    
    liberties = 0
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        new_x, new_y = x + dx, y + dy
        liberties += count_liberties(board, new_x, new_y, visited)
    
    return liberties

def find_best_capture():
    board = create_board()
    best_move = None
    max_capture = 0
    
    # Try each empty position
    for i in range(12):
        for j in range(12):
            if board[i][j] == '.':
                # Simulate placing a black stone
                board[i][j] = 'X'
                # Check surrounding white groups
                for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
                    x, y = i + dx, j + dy
                    if 0 <= x < 12 and 0 <= y < 12 and board[x][y] == 'O':
                        liberties = count_liberties(board, x, y)
                        if liberties == 0:  # This move would capture
                            stones_captured = sum(1 for a in range(12) for b in range(12) 
                                               if board[a][b] == 'O')
                            if stones_captured > max_capture:
                                max_capture = stones_captured
                                best_move = (j+1, 12-i)  # Convert to game coordinates
                board[i][j] = '.'  # Reset the position
    
    if best_move:
        col = chr(ord('A') + best_move[0]-1)
        row = str(best_move[1])
        return f"{col}{row}", max_capture

move, captures = find_best_capture()
print(f"Best move: {move} (captures {captures} stones)")