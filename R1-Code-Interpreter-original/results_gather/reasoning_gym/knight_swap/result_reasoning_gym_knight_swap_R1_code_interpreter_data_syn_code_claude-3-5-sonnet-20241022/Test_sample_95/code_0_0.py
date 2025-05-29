from collections import deque

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 3]

def state_to_string(black_pos, white_pos, empty_pos):
    return ','.join([f"B:{','.join(map(str,black_pos))}", 
                    f"W:{','.join(map(str,white_pos))}", 
                    f"E:{','.join(map(str,empty_pos))}"])

def is_target_state(black_pos, white_pos):
    target_black = {(3,0), (3,1)}  # D1, D2
    target_white = {(1,0), (1,2)}  # B1, B3
    return set(black_pos) == target_black and set(white_pos) == target_white

def solve_knight_swap():
    # Initial positions
    black_pos = [(1,0), (1,2)]  # B1, B3
    white_pos = [(3,0), (3,1)]  # D1, D2
    empty_pos = [(0,0), (0,2), (2,0), (2,2), (1,1)]  # Empty squares
    
    queue = deque([(black_pos, white_pos, empty_pos, [], True)])  # True = Black's turn
    visited = {state_to_string(black_pos, white_pos, empty_pos)}
    
    while queue:
        black_pos, white_pos, empty_pos, moves, black_turn = queue.popleft()
        
        if is_target_state(black_pos, white_pos):
            return moves
            
        current_pieces = black_pos if black_turn else white_pos
        for i, piece in enumerate(current_pieces):
            for new_pos in get_knight_moves(piece):
                if new_pos in empty_pos:
                    # Make move
                    new_black = list(black_pos)
                    new_white = list(white_pos)
                    new_empty = list(empty_pos)
                    
                    if black_turn:
                        new_empty[new_empty.index(new_pos)] = black_pos[i]
                        new_black[i] = new_pos
                    else:
                        new_empty[new_empty.index(new_pos)] = white_pos[i]
                        new_white[i] = new_pos
                    
                    state = state_to_string(new_black, new_white, new_empty)
                    if state not in visited:
                        visited.add(state)
                        move = f"{'B' if black_turn else 'w'},{chr(piece[0]+65)}{piece[1]+1},{chr(new_pos[0]+65)}{new_pos[1]+1}"
                        queue.append((new_black, new_white, new_empty, moves + [move], not black_turn))
    
    return None

solution = solve_knight_swap()
print(solution if solution else "No")