def print_board():
    board = [
        ['O', '.', '.', '.', '.', '.', '.', '.', '.', 'X'],  # 10
        ['.', '.', '.', 'O', '.', '.', '.', '.', '.', '.'],  # 9
        ['.', 'O', '.', '.', 'O', '.', '.', '.', '.', '.'],  # 8
        ['X', 'O', 'X', '.', '.', '.', '.', '.', '.', '.'],  # 7
        ['O', 'O', '.', 'O', '.', 'X', '.', '.', '.', '.'],  # 6
        ['X', 'O', 'X', 'O', '.', '.', '.', '.', '.', '.'],  # 5
        ['.', 'X', '.', '.', '.', 'O', '.', '.', '.', '.'],  # 4
        ['.', '.', '.', '.', '.', '.', '.', 'O', '.', '.'],  # 3
        ['.', '.', '.', 'X', '.', '.', '.', '.', 'X', '.'],  # 2
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 1
    ]
    
    # Check liberties of white group
    def count_liberties(x, y, checked=None):
        if checked is None:
            checked = set()
        
        if x < 0 or x >= 10 or y < 0 or y >= 10:
            return 0
        
        if board[y][x] == '.':
            return 1
        
        if board[y][x] != 'O' or (x,y) in checked:
            return 0
            
        checked.add((x,y))
        return (count_liberties(x+1, y, checked) +
                count_liberties(x-1, y, checked) +
                count_liberties(x, y+1, checked) +
                count_liberties(x, y-1, checked))

    # Check key points
    key_points = [(2,6), (2,5), (1,6)]  # C6, C5, B6
    for x, y in key_points:
        if board[y][x] == '.':
            print(f"Testing point {chr(65+x)}{10-y}")
            # Make move
            board[y][x] = 'X'
            liberties = count_liberties(1, 4)  # Check B5 group
            print(f"Liberties after move: {liberties}")
            # Undo move
            board[y][x] = '.'

print_board()