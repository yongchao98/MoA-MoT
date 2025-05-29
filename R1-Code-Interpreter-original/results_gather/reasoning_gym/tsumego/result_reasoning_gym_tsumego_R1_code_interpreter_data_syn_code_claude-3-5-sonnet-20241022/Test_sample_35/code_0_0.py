def create_board():
    board = {}
    # Empty board with . representing empty points
    for row in range(1, 10):
        for col in 'ABCDEFGHI':
            board[col + str(row)] = '.'
    
    # Add White stones (O)
    white_stones = ['G9', 'G8', 'C8', 'C7', 'B6', 'C6', 'D6', 'A6', 'C5', 'E3', 'H3', 'D2', 'A1']
    for pos in white_stones:
        board[pos] = 'O'
    
    # Add Black stones (X)
    black_stones = ['D7', 'B6', 'E6', 'A5', 'B5', 'C4', 'E2']
    for pos in black_stones:
        board[pos] = 'X'
    
    return board

def get_adjacent(pos):
    col, row = pos[0], int(pos[1])
    adj = []
    # Check all four directions
    if col != 'A': adj.append(chr(ord(col)-1) + str(row))  # Left
    if col != 'I': adj.append(chr(ord(col)+1) + str(row))  # Right
    if row > 1: adj.append(col + str(row-1))  # Down
    if row < 9: adj.append(col + str(row+1))  # Up
    return adj

def find_group_liberties(board, start_pos):
    if board[start_pos] == '.':
        return set(), set()
    
    color = board[start_pos]
    group = set([start_pos])
    liberties = set()
    visited = set([start_pos])
    to_check = [start_pos]
    
    while to_check:
        pos = to_check.pop()
        for adj in get_adjacent(pos):
            if adj not in visited:
                visited.add(adj)
                if board[adj] == color:
                    group.add(adj)
                    to_check.append(adj)
                elif board[adj] == '.':
                    liberties.add(adj)
    
    return group, liberties

def find_capturable_moves(board):
    potential_moves = set()
    
    # Check each position on the board
    for pos in board:
        if board[pos] == 'O':  # For each White stone
            group, liberties = find_group_liberties(board, pos)
            if len(liberties) == 1:  # If group has only one liberty
                potential_moves.update(liberties)
    
    return potential_moves

# Create and analyze the board
board = create_board()
capturable_moves = find_capturable_moves(board)
print("Possible capturing moves:", sorted(list(capturable_moves)))