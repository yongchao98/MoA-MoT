def print_board_with_liberties():
    board = [
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', 'O', '.', 'O', 'O', 'O', '.'],
        ['.', '.', '.', 'X', '.', 'X', '.', '.', '.', '.'],
        ['.', '.', 'X', 'O', 'O', 'O', 'X', '.', '.', '.'],
        ['.', '.', '.', 'X', 'O', 'X', '.', '.', '.', '.'],
        ['.', '.', '.', 'X', 'X', '.', '.', '.', '.', '.'],
        ['.', '.', '.', 'O', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', 'O', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
        ['O', 'O', '.', '.', '.', '.', '.', 'X', '.', '.']
    ]
    
    # Mark liberties of the central white group with '?'
    liberties = set()
    white_group = [(4,4), (3,4), (3,3), (3,5)]  # E6, E7, D7, F7
    
    for x, y in white_group:
        for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
            nx, ny = x + dx, y + dy
            if 0 <= nx < 10 and 0 <= ny < 10:
                if board[nx][ny] == '.':
                    liberties.add((nx,ny))
    
    # Print board with liberties
    print("   A B C D E F G H I J")
    for i in range(10):
        row = f"{10-i:2d}"
        for j in range(10):
            if (i,j) in liberties:
                row += " ?"
            else:
                row += " " + board[i][j]
        print(row)
    
    print("\nPotential liberty points:", sorted([(chr(y+65), 10-x) for x,y in liberties]))

print_board_with_liberties()