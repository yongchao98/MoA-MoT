def print_board_analysis(board):
    # Convert board representation for easier visualization
    symbols = {'.': '⋅', 'X': '●', 'O': '○'}
    liberties = {}
    
    def count_liberties(x, y, stone):
        if (x, y) in liberties:
            return liberties[(x, y)]
        
        visited = set()
        group_liberties = set()
        
        def explore(cx, cy):
            if (cx, cy) in visited:
                return
            if not (0 <= cx < 10 and 0 <= cy < 10):
                return
            if board[cy][cx] == '.':
                group_liberties.add((cx, cy))
                return
            if board[cy][cx] != stone:
                return
                
            visited.add((cx, cy))
            for nx, ny in [(cx+1,cy), (cx-1,cy), (cx,cy+1), (cx,cy-1)]:
                explore(nx, ny)
        
        explore(x, y)
        liberties[(x, y)] = len(group_liberties)
        return len(group_liberties)
    
    # Print board with liberty count for white groups
    for y in range(10):
        row = f"{10-y:2d} "
        for x in range(10):
            if board[y][x] == 'O':
                libs = count_liberties(x, y, 'O')
                row += f"{symbols[board[y][x]]}({libs}) "
            else:
                row += f"{symbols[board[y][x]]}    "
        print(row)
    print("   A    B    C    D    E    F    G    H    I    J")
    
    # Analyze potential capturing moves
    for y in range(10):
        for x in range(10):
            if board[y][x] == 'O':
                libs = count_liberties(x, y, 'O')
                if libs <= 2:  # Focus on stones with 2 or fewer liberties
                    print(f"White group at {chr(65+x)}{10-y} has {libs} liberties")

board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', 'O', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'O', '.', '.', '.', 'X', '.', '.'],
    ['.', '.', '.', '.', 'X', '.', '.', '.', 'X', 'O'],
    ['.', '.', '.', '.', '.', '.', 'X', 'O', 'X', 'X'],
    ['O', '.', '.', '.', '.', '.', '.', 'X', 'O', 'O'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', 'O', '.']
]

print_board_analysis(board)