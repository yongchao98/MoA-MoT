from collections import deque

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 3]

def state_to_string(whites, blacks, turn):
    return ','.join([f'{x},{y}' for x,y in sorted(whites)]) + '|' + \
           ','.join([f'{x},{y}' for x,y in sorted(blacks)]) + f'|{turn}'

def string_to_pos(s):
    x, y = map(int, s.split(','))
    return (x, y)

def solve_knight_swap():
    # Initial state
    init_whites = {(2,2), (0,1)}  # C3, A2
    init_blacks = {(0,0), (2,0)}  # A1, C1
    target_blacks = init_whites
    target_whites = init_blacks
    
    queue = deque([(init_whites, init_blacks, [], 'w')])
    visited = {state_to_string(init_whites, init_blacks, 'w')}
    
    while queue:
        whites, blacks, moves, turn = queue.pop()
        
        # Check if reached goal
        if whites == target_whites and blacks == target_blacks:
            return moves
        
        # Generate next moves
        current_pieces = whites if turn == 'w' else blacks
        other_pieces = blacks if turn == 'w' else whites
        next_turn = 'b' if turn == 'w' else 'w'
        
        for piece in current_pieces:
            for new_pos in get_knight_moves(piece):
                if new_pos not in whites and new_pos not in blacks:
                    new_current = (current_pieces - {piece}) | {new_pos}
                    new_whites = new_current if turn == 'w' else whites
                    new_blacks = new_current if turn == 'b' else blacks
                    
                    state = state_to_string(new_whites, new_blacks, next_turn)
                    if state not in visited:
                        visited.add(state)
                        new_moves = moves + [f"{turn},{chr(piece[0]+65)}{piece[1]+1},{chr(new_pos[0]+65)}{new_pos[1]+1}"]
                        queue.appendleft((new_whites, new_blacks, new_moves, next_turn))
    
    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")