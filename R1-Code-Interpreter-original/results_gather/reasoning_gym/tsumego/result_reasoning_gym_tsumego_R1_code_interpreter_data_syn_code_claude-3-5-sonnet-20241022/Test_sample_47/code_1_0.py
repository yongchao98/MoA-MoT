def check_capture_size(board, x, y, visited=None):
    if visited is None:
        visited = set()
    
    if not (0 <= x < 10 and 0 <= y < 10) or board[x][y] != 'O':
        return 0
    
    if (x, y) in visited:
        return 0
        
    visited.add((x, y))
    count = 1
    
    # Count connected stones
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        new_x, new_y = x + dx, y + dy
        count += check_capture_size(board, new_x, new_y, visited)
    
    return count

def analyze_moves(board):
    moves = {}
    for i in range(10):
        for j in range(10):
            if board[i][j] == '.':
                # Try each empty point
                capture_potential = 0
                for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
                    new_x, new_y = i + dx, j + dy
                    if 0 <= new_x < 10 and 0 <= new_y < 10 and board[new_x][new_y] == 'O':
                        # Check if this move would capture stones
                        liberties = get_liberties(board, new_x, new_y)
                        if liberties == 1:  # Can be captured
                            capture_potential += check_capture_size(board, new_x, new_y)
                
                if capture_potential > 0:
                    move = (chr(ord('A') + j), 10 - i)
                    moves[move] = capture_potential

    return moves

board = create_board()
capturing_moves = analyze_moves(board)
print("Moves and number of stones they can capture:", capturing_moves)