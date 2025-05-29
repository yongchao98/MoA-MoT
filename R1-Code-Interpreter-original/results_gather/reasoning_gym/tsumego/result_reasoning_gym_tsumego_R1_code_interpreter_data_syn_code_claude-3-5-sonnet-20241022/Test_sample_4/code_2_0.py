def create_board():
    board = {}
    for i in range(1, 12):
        for c in 'ABCDEFGHIJK':
            board[f"{c}{i}"] = '.'
    
    # Set white stones (O)
    white_stones = ['A11', 'H9', 'G6', 'H5', 'I5', 'D4', 'H4', 'I4', 'I3']
    for pos in white_stones:
        board[pos] = 'O'
    
    # Set black stones (X)
    black_stones = ['B7', 'I7', 'I6', 'C5', 'J5', 'G4', 'J4', 'H3', 'J3', 'I2', 'C1', 'G1']
    for pos in black_stones:
        board[pos] = 'X'
    
    return board

def check_group_capture(board, move):
    # Simulate the move
    test_board = board.copy()
    test_board[move] = 'X'
    
    # Check if any adjacent white stones would be captured
    def get_group_liberties(pos, checked=None):
        if checked is None:
            checked = set()
        if pos in checked:
            return set()
        checked.add(pos)
        liberties = set()
        for adj in [f"{chr(ord(pos[0])+dx)}{int(pos[1])+dy}" 
                   for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]
                   if chr(ord(pos[0])+dx) in 'ABCDEFGHIJK' and 
                   1 <= int(pos[1])+dy <= 11]:
            if adj not in test_board or test_board[adj] == '.':
                liberties.add(adj)
            elif test_board[adj] == 'O':
                liberties.update(get_group_liberties(adj, checked))
        return liberties

    # Check adjacent positions for white stones
    for adj in [f"{chr(ord(move[0])+dx)}{int(move[1])+dy}" 
                for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]
                if chr(ord(move[0])+dx) in 'ABCDEFGHIJK' and 
                1 <= int(move[1])+dy <= 11]:
        if adj in test_board and test_board[adj] == 'O':
            if not get_group_liberties(adj):
                return True
    return False

board = create_board()
move = 'G5'
will_capture = check_group_capture(board, move)
print(f"Playing at {move} will capture white stones: {will_capture}")