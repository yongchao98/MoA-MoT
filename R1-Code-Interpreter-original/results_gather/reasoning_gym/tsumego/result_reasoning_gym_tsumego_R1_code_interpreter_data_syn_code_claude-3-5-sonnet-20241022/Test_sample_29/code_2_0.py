def check_sequence():
    # Initialize board state
    board = [['.'] * 12 for _ in range(12)]
    
    # Initial position
    black = [(9,9), (9,10), (9,11), (9,12),  # J10,K10,L10,M10
            (7,9), (6,9), (6,10)]           # J8,J7,K7
    white = [(8,9), (8,10), (8,11), (7,10)]  # J9,K9,L9,K8
    
    # Place stones
    for x, y in black:
        board[x][y] = 'X'
    for x, y in white:
        board[x][y] = 'O'
    
    # Try L8 move
    board[7][11] = 'X'  # L8 move
    
    # Count liberties of white group after move
    def count_group_liberties(x, y, visited=None):
        if visited is None:
            visited = set()
        
        if (x, y) not in visited and 0 <= x < 12 and 0 <= y < 12:
            if board[x][y] == '.':
                return 1
            if board[x][y] == 'O':
                visited.add((x, y))
                liberties = 0
                for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
                    liberties += count_group_liberties(x+dx, y+dy, visited)
                return liberties
        return 0
    
    # Check liberties of white group at K8
    liberties = count_group_liberties(7,10)
    print(f"White group liberties after Black plays L8: {liberties}")
    
    # Check if white can escape
    escape_points = []
    for x, y in [(7,11), (8,11), (6,10)]:  # Adjacent points
        if 0 <= x < 12 and 0 <= y < 12 and board[x][y] == '.':
            escape_points.append((chr(y + ord('A')), x + 1))
    
    print(f"Possible escape points for White: {escape_points}")
    return liberties

liberties = check_sequence()