def count_liberties(board, x, y):
    if board[y][x] == '.':
        return -1
    
    color = board[y][x]
    visited = set()
    liberties = set()
    
    def flood_fill(x, y):
        if (x, y) in visited:
            return
        if x < 0 or x >= 10 or y < 0 or y >= 10:
            return
        if board[y][x] != color and board[y][x] != '.':
            return
        
        if board[y][x] == '.':
            liberties.add((x, y))
            return
            
        visited.add((x, y))
        flood_fill(x+1, y)
        flood_fill(x-1, y)
        flood_fill(x, y+1)
        flood_fill(x, y-1)
    
    flood_fill(x, y)
    return len(liberties)

# Create the board
board = [
    ['.', '.', 'O', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'O', '.'],
    ['.', '.', 'X', '.', '.', '.', 'O', '.', 'X', '.'],
    ['.', '.', 'X', '.', '.', '.', 'X', 'X', '.', '.'],
    ['.', '.', '.', '.', 'O', 'X', 'O', 'O', 'X', '.'],
    ['.', '.', '.', '.', '.', 'X', 'O', '.', '.', '.'],
    ['.', 'X', '.', '.', '.', '.', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'O', '.', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', '.']
]

# Check liberties for white groups
critical_positions = []
for y in range(10):
    for x in range(10):
        if board[y][x] == 'O':
            liberties = count_liberties(board, x, y)
            if liberties <= 2:  # Critical if 2 or fewer liberties
                critical_positions.append((chr(ord('A') + x), 10-y, liberties))

print("Critical white positions (position, liberties):")
for pos in critical_positions:
    print(f"{pos[0]}{pos[1]}: {pos[2]} liberties")