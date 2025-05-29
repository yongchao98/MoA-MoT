def print_possible_moves(board):
    def get_knight_moves(pos):
        x, y = pos
        moves = [
            (x+2, y+1), (x+2, y-1), (x-2, y+1), (x-2, y-1),
            (x+1, y+2), (x+1, y-2), (x-1, y+2), (x-1, y-2)
        ]
        return [(x, y) for x, y in moves if 0 <= x < 4 and 0 <= y < 3]
    
    # Initial board
    board = [
        ['.', '.', 'B', 'w'],
        ['.', 'B', '.', '.'],
        ['.', 'w', '.', '.']
    ]
    
    # Find white knights (first to move)
    white_moves = []
    for i in range(4):
        for j in range(3):
            if board[j][i] == 'w':
                pos = (i, j)
                valid_moves = [(x,y) for x,y in get_knight_moves(pos) if board[y][x] == '.']
                if valid_moves:
                    white_moves.append((f"{chr(65+pos[0])}{pos[1]+1}", 
                                     [f"{chr(65+x)}{y+1}" for x,y in valid_moves]))
    
    print("Possible first moves for white knights:")
    for start, moves in white_moves:
        print(f"From {start}: {moves}")

print_possible_moves()