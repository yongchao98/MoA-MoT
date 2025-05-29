from collections import deque

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 3]

def board_to_string(knights):
    return ','.join(sorted([f"{color}{x}{y}" for color, (x,y) in knights.items()]))

def is_target_reached(knights):
    white_positions = set(pos for color, pos in knights.items() if color.startswith('w'))
    black_original = {(1,2), (2,0)}  # B3, C1
    return white_positions == black_original

def solve_knight_swap():
    initial_knights = {
        'w1': (2,1), 'w2': (3,1),  # C2, D2
        'b1': (1,2), 'b2': (2,0)   # B3, C1
    }
    
    visited = set()
    queue = deque([(initial_knights, [], True)])
    min_solution = None
    
    while queue:
        knights, moves, is_white_turn = queue.popleft()
        state_key = board_to_string(knights)
        
        if state_key in visited:
            continue
        visited.add(state_key)
        
        if is_target_reached(knights) and len(moves) > 0:
            if min_solution is None or len(moves) < len(min_solution):
                min_solution = moves
            continue
            
        if min_solution and len(moves) >= len(min_solution):
            continue
            
        current_pieces = [k for k in knights.keys() if k.startswith('w' if is_white_turn else 'b')]
        
        for piece in current_pieces:
            x, y = knights[piece]
            for new_x, new_y in get_knight_moves((x, y)):
                if not any((px,py) == (new_x,new_y) for px,py in knights.values()):
                    new_knights = knights.copy()
                    new_knights[piece] = (new_x, new_y)
                    move = f"{'w' if is_white_turn else 'b'},{chr(65+x)}{y+1},{chr(65+new_x)}{new_y+1}"
                    queue.append((new_knights, moves + [move], not is_white_turn))
    
    return min_solution

# Solve and print result
result = solve_knight_swap()
if result:
    print(result)
else:
    print("No")