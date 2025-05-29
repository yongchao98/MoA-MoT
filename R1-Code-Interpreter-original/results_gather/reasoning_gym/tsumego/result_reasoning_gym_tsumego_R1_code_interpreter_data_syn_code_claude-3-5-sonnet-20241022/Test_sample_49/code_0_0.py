def create_board():
    board = {}
    # Convert the given position to a dictionary format
    position = """
   . . . . . . . . . . . . .
   . . . . . . . . . . . . .
   . . X . . . . . . . . . .
   . . . . O . . . . . . . .
   . . . . . . . . . . . . .
   . . . . O . . . . . . . .
   . . . . X . . O . . . . X
   . . . X . O . . . . . . .
   . . O O O O X . . . . . .
   . X . X O X . . . . . . .
   . . . . X . . . . . . . .
   . . . . . . . . . . . . .
   . . . . O . . . . X . X .
    """
    rows = position.strip().split('\n')
    for i, row in enumerate(rows):
        stones = row.strip().split()
        for j, stone in enumerate(stones):
            if stone in ['O', 'X']:
                board[(j, 12-i)] = stone
    return board

def get_liberties(board, pos):
    x, y = pos
    liberties = []
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        new_pos = (x+dx, y+dy)
        if 0 <= new_pos[0] < 13 and 0 <= new_pos[1] < 13:
            if new_pos not in board:
                liberties.append(new_pos)
    return liberties

def find_white_groups(board):
    white_stones = [(x,y) for (x,y) in board if board[(x,y)] == 'O']
    groups = []
    visited = set()
    
    def find_group(pos):
        group = set([pos])
        stack = [pos]
        while stack:
            current = stack.pop()
            for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
                new_pos = (current[0]+dx, current[1]+dy)
                if (new_pos in board and 
                    board[new_pos] == 'O' and 
                    new_pos not in group):
                    group.add(new_pos)
                    stack.append(new_pos)
        return group

    for stone in white_stones:
        if stone not in visited:
            group = find_group(stone)
            visited.update(group)
            groups.append(group)
    
    return groups

def analyze_position():
    board = create_board()
    white_groups = find_white_groups(board)
    
    # Analyze each white group's liberties
    critical_moves = []
    for group in white_groups:
        group_liberties = set()
        for stone in group:
            group_liberties.update(get_liberties(board, stone))
        if len(group_liberties) <= 2:
            critical_moves.append((group_liberties, len(group)))
    
    # Find the move that could capture the most stones
    best_move = None
    max_capture = 0
    for liberties, group_size in critical_moves:
        if len(liberties) == 1:
            liberty = list(liberties)[0]
            if group_size > max_capture:
                max_capture = group_size
                best_move = liberty
    
    if best_move:
        col = chr(ord('A') + best_move[0])
        row = best_move[1] + 1
        print(f"Best move: {col}{row} can capture {max_capture} stones")

analyze_position()