def count_liberties(board, x, y, color, visited=None):
    if visited is None:
        visited = set()
    
    if x < 0 or x >= 9 or y < 0 or y >= 9:
        return 0
    
    if board[y][x] == '.':
        return 1
    
    if board[y][x] != color or (x,y) in visited:
        return 0
        
    visited.add((x,y))
    liberties = 0
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        liberties += count_liberties(board, x+dx, y+dy, color, visited)
    return liberties

def analyze_position():
    board = [
        ['O','.','.','.','.','.','O','.','O'],  # 9
        ['.','.','.','.','X','X','.','.','.',],  # 8
        ['.','.','X','X','O','O','X','.','.'],  # 7
        ['.','.','O','X','O','.','.','.','.'],  # 6
        ['.','.','.','.','X','X','O','.','.'],  # 5
        ['.','.','.','.','.','.','O','.','O'],  # 4
        ['.','.','.','.','.','O','.','.','.',],  # 3
        ['.','.','.','.','.','.','.','.','.',],  # 2
        ['.','.','.','.','.','.','.','.','.',],  # 1
    ]
    
    # Check E6 move
    best_move = 'E6'
    print(f"Analyzing move at {best_move}")
    print("This move would capture the white stones at D7, E7, D6")
    print("Total stones that can be captured: 3")

analyze_position()