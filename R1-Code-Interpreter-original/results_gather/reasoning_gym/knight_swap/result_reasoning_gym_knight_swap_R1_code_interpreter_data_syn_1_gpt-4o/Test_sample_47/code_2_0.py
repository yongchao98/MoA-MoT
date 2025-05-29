from collections import deque

# Define the initial positions of the knights
initial_white_knights = {'A2': 'w', 'D1': 'w'}
initial_black_knights = {'A1': 'B', 'B2': 'B'}

# Define the target positions
target_white = {'A1', 'B2'}
target_black = {'A2', 'D1'}

# Function to get possible moves for a knight at a given position
def get_knight_moves(pos):
    col, row = pos
    col_idx = ord(col) - ord('A')
    row_idx = int(row) - 1
    potential_moves = [
        (col_idx + 2, row_idx + 1), (col_idx + 2, row_idx - 1),
        (col_idx - 2, row_idx + 1), (col_idx - 2, row_idx - 1),
        (col_idx + 1, row_idx + 2), (col_idx + 1, row_idx - 2),
        (col_idx - 1, row_idx + 2), (col_idx - 1, row_idx - 2)
    ]
    valid_moves = []
    for c, r in potential_moves:
        if 0 <= c < 4 and 0 <= r < 3:
            valid_moves.append(chr(c + ord('A')) + str(r + 1))
    return valid_moves

# Function to check if the current configuration matches the target
def is_target_reached(white_knights, black_knights):
    return set(white_knights.keys()) == target_white and set(black_knights.keys()) == target_black

# BFS to find the sequence of moves
def bfs_knight_swap():
    queue = deque()
    queue.append((initial_white_knights, initial_black_knights, [], 'w'))
    visited = set()

    while queue:
        white_knights, black_knights, moves, current_turn = queue.popleft()
        state = (frozenset(white_knights.items()), frozenset(black_knights.items()), current_turn)
        
        if state in visited:
            continue
        visited.add(state)

        if is_target_reached(white_knights, black_knights):
            return moves

        if current_turn == 'w':
            for pos in list(white_knights.keys()):
                for move in get_knight_moves(pos):
                    if move not in white_knights and move not in black_knights:
                        new_white_knights = white_knights.copy()
                        new_white_knights[move] = new_white_knights.pop(pos)
                        queue.append((new_white_knights, black_knights, moves + [f"w,{pos},{move}"], 'B'))
        else:
            for pos in list(black_knights.keys()):
                for move in get_knight_moves(pos):
                    if move not in white_knights and move not in black_knights:
                        new_black_knights = black_knights.copy()
                        new_black_knights[move] = new_black_knights.pop(pos)
                        queue.append((white_knights, new_black_knights, moves + [f"B,{pos},{move}"], 'w'))

    return "No"

# Output the result
result = bfs_knight_swap()
print(result)