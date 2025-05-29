def create_board():
    board = {}
    # Initialize empty board
    for row in range(1, 14):
        for col in 'ABCDEFGHIJKLM':
            board[f"{col}{row}"] = '.'
    
    # Add Black stones (X)
    black_stones = ["F10", "H10", "H9", "B5", "G4", "C3", "E3", "B2", "F2", "I2", "C1", "E1"]
    for pos in black_stones:
        board[pos] = 'X'
    
    # Add White stones (O)
    white_stones = ["I13", "I9", "K9", "D4", "B2", "C2", "D2", "E2", "D1", "H1"]
    for pos in white_stones:
        board[pos] = 'O'
    
    return board

def get_adjacent(pos):
    col, row = pos[0], int(pos[1:])
    adj = []
    # Check all four directions
    if col > 'A': adj.append(f"{chr(ord(col)-1)}{row}")
    if col < 'M': adj.append(f"{chr(ord(col)+1)}{row}")
    if row > 1: adj.append(f"{col}{row-1}")
    if row < 13: adj.append(f"{col}{row+1}")
    return adj

def find_group_liberties(board, start_pos):
    if board[start_pos] == '.': return set(), set()
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

def analyze_position():
    board = create_board()
    white_groups = []
    
    # Find all white groups and their liberties
    checked = set()
    for pos in board:
        if board[pos] == 'O' and pos not in checked:
            group, liberties = find_group_liberties(board, pos)
            checked.update(group)
            white_groups.append((group, liberties))
    
    # Print analysis
    for i, (group, liberties) in enumerate(white_groups):
        print(f"White group {i+1}: stones={sorted(group)}, liberties={sorted(liberties)}")

analyze_position()