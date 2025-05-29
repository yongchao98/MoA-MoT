def count_liberties(board, x, y, visited=None):
    if visited is None:
        visited = set()
    
    if x < 0 or x >= 9 or y < 0 or y >= 9:
        return 0
    
    if board[y][x] == '.':
        return 1
    
    if (x, y) in visited:
        return 0
        
    if board[y][x] == 'O':
        visited.add((x, y))
        liberties = 0
        # Check all four directions
        for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
            new_x, new_y = x + dx, y + dy
            if 0 <= new_x < 9 and 0 <= new_y < 9:
                if board[new_y][new_x] == '.':
                    liberties += 1
                elif board[new_y][new_x] == 'O':
                    liberties += count_liberties(board, new_x, new_y, visited)
        return liberties
    return 0

# Create the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', '.', '.', 'O', '.', '.'],
    ['.', 'O', 'X', 'X', '.', '.', 'O', '.', '.'],
    ['.', 'X', 'O', 'O', 'X', 'O', '.', 'O', '.'],
    ['.', 'X', 'O', '.', '.', '.', 'O', '.', '.'],
    ['O', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', 'X', '.', '.', '.', 'X'],
    ['.', 'O', '.', '.', '.', '.', '.', '.', '.']
]

# Check each empty point for potential capturing moves
for y in range(9):
    for x in range(9):
        if board[y][x] == '.':
            # Try placing a black stone here
            board[y][x] = 'X'
            # Check surrounding white groups
            for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
                new_x, new_y = x + dx, y + dy
                if 0 <= new_x < 9 and 0 <= new_y < 9:
                    if board[new_y][new_x] == 'O':
                        liberties = count_liberties(board, new_x, new_y)
                        if liberties == 0:
                            print(f"Capturing move found at {chr(65+x)}{9-y}")
            board[y][x] = '.'  # Reset the position