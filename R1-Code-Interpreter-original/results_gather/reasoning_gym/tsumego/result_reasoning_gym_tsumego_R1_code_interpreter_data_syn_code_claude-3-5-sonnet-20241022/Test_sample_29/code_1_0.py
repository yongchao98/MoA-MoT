def count_liberties(x, y, board, checked):
    if (x, y) in checked or x < 0 or y < 0 or x >= 12 or y >= 12:
        return 0
    if board[x][y] == '.':
        return 1
    return 0

def analyze_position():
    board = [['.'] * 12 for _ in range(12)]
    # Place the stones accurately
    white_positions = [(8,9), (8,10), (8,11), (7,10)]  # J9,K9,L9,K8
    black_positions = [(9,9), (9,10), (9,11), (9,12),  # J10,K10,L10,M10
                      (7,9), (6,9), (6,10)]            # J8,J7,K7
    
    # Count liberties for white group
    total_liberties = 0
    checked = set()
    
    for x, y in white_positions:
        # Check adjacent points
        for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
            total_liberties += count_liberties(x+dx, y+dy, board, checked)
            checked.add((x+dx, y+dy))
    
    # Analyze key points
    key_moves = ['L8', 'K7', 'J8']
    print("Potential key moves:", key_moves)
    print("White group liberties:", total_liberties)
    
    return key_moves[0]  # L8 is the vital point

result = analyze_position()
print("\nCritical move:", result)