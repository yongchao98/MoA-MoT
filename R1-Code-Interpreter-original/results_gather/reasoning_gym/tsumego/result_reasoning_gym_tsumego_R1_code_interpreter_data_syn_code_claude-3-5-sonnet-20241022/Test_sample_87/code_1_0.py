def check_group_capture(x, y, board):
    moves = [(x+1,y), (x-1,y), (x,y+1), (x,y-1)]
    white_group = set()
    liberties = set()
    
    def find_group(x, y, stone_type):
        if x < 0 or x >= 11 or y < 0 or y >= 11:
            return
        if board[y][x] != stone_type:
            if board[y][x] == '.':
                liberties.add((x,y))
            return
        if (x,y) in white_group:
            return
        white_group.add((x,y))
        for dx, dy in [(1,0), (-1,0), (0,1), (0,-1)]:
            find_group(x+dx, y+dy, stone_type)
    
    find_group(x, y, 'O')
    
    print(f"White group from {chr(65+x)}{y+1}:")
    print(f"Group size: {len(white_group)}")
    print("Liberties at:", end=" ")
    for lx, ly in liberties:
        print(f"{chr(65+lx)}{ly+1}", end=" ")
    print("\n")
    return liberties

def analyze_critical_position():
    board = [
        ['.', '.', '.', '.', '.', '.', '.', 'X', '.', '.', '.'],  # 1
        ['.', '.', '.', '.', '.', '.', 'X', 'O', '.', '.', '.'],  # 2
        ['.', '.', '.', '.', '.', '.', 'X', 'O', 'O', 'X', '.'],  # 3
        ['.', '.', '.', '.', '.', '.', 'O', 'X', 'X', '.', '.'],  # 4
        ['O', '.', '.', '.', '.', '.', '.', 'O', '.', '.', 'O'],  # 5
        ['.', '.', '.', 'X', '.', '.', '.', '.', '.', '.', '.'],  # 6
        ['.', '.', '.', 'X', '.', '.', '.', '.', '.', '.', '.'],  # 7
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 8
        ['.', '.', '.', '.', 'O', '.', '.', '.', 'O', '.', '.'],  # 9
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 10
        ['.', '.', '.', 'X', '.', '.', 'O', '.', 'O', '.', '.'],  # 11
    ]
    
    # Analyze the white group starting from H3
    liberties = check_group_capture(7, 2, board)  # H3
    return liberties

analyze_critical_position()