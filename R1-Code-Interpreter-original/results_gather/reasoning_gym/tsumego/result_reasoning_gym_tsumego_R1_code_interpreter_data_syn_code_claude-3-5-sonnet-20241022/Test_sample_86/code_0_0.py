def create_board():
    board = {}
    for i in range(1, 11):
        for j in 'ABCDEFGHIJ':
            board[j + str(i)] = '.'
    
    # Set black stones (X)
    black_stones = ['A6', 'H4', 'G3', 'I3', 'F2', 'J2', 'G1', 'I1']
    for pos in black_stones:
        board[pos] = 'X'
    
    # Set white stones (O)
    white_stones = ['I10', 'J9', 'B6', 'F5', 'E4', 'H3', 'A2', 'C2', 'G2', 'H2', 'C1', 'F1', 'H1']
    for pos in white_stones:
        board[pos] = 'O'
    
    return board

def get_adjacent(pos):
    col, row = pos[0], int(pos[1:])
    adj = []
    if col > 'A': adj.append(chr(ord(col)-1) + str(row))
    if col < 'J': adj.append(chr(ord(col)+1) + str(row))
    if row > 1: adj.append(col + str(row-1))
    if row < 10: adj.append(col + str(row+1))
    return adj

def find_group(board, start_pos, stone_type, visited=None):
    if visited is None:
        visited = set()
    
    if start_pos not in board or board[start_pos] != stone_type:
        return set()
    
    group = {start_pos}
    visited.add(start_pos)
    
    for adj in get_adjacent(start_pos):
        if adj not in visited and adj in board and board[adj] == stone_type:
            group.update(find_group(board, adj, stone_type, visited))
    
    return group

def has_liberties(board, group):
    for pos in group:
        for adj in get_adjacent(pos):
            if adj in board and board[adj] == '.':
                return True
    return False

def count_captures(board, move):
    if board[move] != '.':
        return 0
    
    # Make the move
    test_board = board.copy()
    test_board[move] = 'X'
    
    captured = set()
    for adj in get_adjacent(move):
        if adj in test_board and test_board[adj] == 'O':
            group = find_group(test_board, adj, 'O')
            if not has_liberties(test_board, group):
                captured.update(group)
    
    return len(captured)

def find_best_move():
    board = create_board()
    best_move = None
    max_captures = 0
    
    for pos in board:
        if board[pos] == '.':
            captures = count_captures(board, pos)
            if captures > max_captures:
                max_captures = captures
                best_move = pos
    
    print(f"Best move: {best_move}")
    print(f"Captures: {max_captures}")

find_best_move()