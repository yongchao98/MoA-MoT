def create_board():
    board = [['.'] * 10 for _ in range(10)]
    # Add black stones (X)
    black_stones = [(7,1), (7,2), (7,3), (6,3), (5,0), (5,2), (4,0), (4,1), (2,6), (2,8)]
    # Add white stones (O)
    white_stones = [(5,1), (5,2), (5,0), (3,0), (3,4), (1,5), (1,7), (0,3), (0,4), (0,9)]
    
    for x, y in black_stones:
        board[x][y] = 'X'
    for x, y in white_stones:
        board[x][y] = 'O'
    
    return board

def count_liberties(board, x, y, visited=None):
    if visited is None:
        visited = set()
    
    if x < 0 or x >= 10 or y < 0 or y >= 10:
        return 0
    
    if board[x][y] == '.':
        return 1
    
    if (x,y) in visited:
        return 0
        
    visited.add((x,y))
    
    liberties = 0
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        new_x, new_y = x + dx, y + dy
        liberties += count_liberties(board, new_x, new_y, visited)
    
    return liberties

def find_vulnerable_groups(board):
    vulnerable = []
    for i in range(10):
        for j in range(10):
            if board[i][j] == 'O':
                liberties = count_liberties(board, i, j)
                if liberties <= 2:  # Groups with 2 or fewer liberties are vulnerable
                    vulnerable.append((i, j, liberties))
    return vulnerable

board = create_board()
vulnerable = find_vulnerable_groups(board)
print("Vulnerable white groups (row, col, liberties):")
for v in vulnerable:
    print(f"Position: ({9-v[0]},{chr(65+v[1])}) - Liberties: {v[2]}")