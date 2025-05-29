def capture_stones(board, move):
    # Convert board to a 2D list
    board = [list(row) for row in board.split('\n')]
    
    # Place the move
    x, y = move
    board[y][x] = 'X'
    
    # Function to check if a group of stones is captured
    def is_captured(x, y, color):
        visited = set()
        stack = [(x, y)]
        captured = True
        
        while stack:
            cx, cy = stack.pop()
            if (cx, cy) in visited:
                continue
            visited.add((cx, cy))
            
            # Check liberties
            for nx, ny in [(cx-1, cy), (cx+1, cy), (cx, cy-1), (cx, cy+1)]:
                if 0 <= nx < len(board[0]) and 0 <= ny < len(board):
                    if board[ny][nx] == '.':
                        captured = False
                    elif board[ny][nx] == color:
                        stack.append((nx, ny))
        
        return captured, visited
    
    # Check all white stones
    captured_stones = set()
    for y in range(len(board)):
        for x in range(len(board[0])):
            if board[y][x] == 'O':
                captured, stones = is_captured(x, y, 'O')
                if captured:
                    captured_stones.update(stones)
    
    return len(captured_stones)

# Define the board
board = """
...........
...........
........X.O
.O.........
..OOXX.....
...XOOXXO..
O..XO......
....X......
..X..O.....
...X.......
..........X
"""

# Move to test
move = (8, 2)  # J9

# Calculate captured stones
captured = capture_stones(board, move)
print(captured)