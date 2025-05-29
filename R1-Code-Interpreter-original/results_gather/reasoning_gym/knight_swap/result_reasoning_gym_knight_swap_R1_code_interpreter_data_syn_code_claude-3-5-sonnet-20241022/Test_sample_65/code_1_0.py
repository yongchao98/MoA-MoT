from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 3]

def pos_to_notation(pos):
    return chr(65 + pos[0]) + str(pos[1] + 1)

def state_to_string(state):
    return ','.join(''.join(row) for row in state)

def is_target_state(state):
    # Check if white knights are at original black positions and vice versa
    target_whites = {(2,2), (3,1)}  # C3, D2
    target_blacks = {(0,2), (1,0)}  # A3, B1
    
    current_whites = set()
    current_blacks = set()
    for i in range(4):
        for j in range(3):
            if state[j][i] == 'w':
                current_whites.add((i,j))
            elif state[j][i] == 'B':
                current_blacks.add((i,j))
    
    return current_whites == target_whites and current_blacks == target_blacks

def solve_knight_swap():
    # Corrected initial state
    initial_state = [
        ['w', '.', 'B', '.'],
        ['.', '.', '.', 'B'],
        ['.', 'w', '.', '.']
    ]
    
    queue = deque([(initial_state, [], True)])  # state, moves, is_black_turn
    seen_states = {state_to_string(initial_state)}
    
    while queue:
        current_state, moves, is_black_turn = queue.popleft()
        
        if is_target_state(current_state):
            return moves
        
        # Find all pieces that can move
        piece = 'B' if is_black_turn else 'w'
        for i in range(4):
            for j in range(3):
                if current_state[2-j][i] == piece:  # Corrected indexing
                    # Try all possible moves for this piece
                    for new_x, new_y in get_knight_moves((i,j)):
                        if current_state[2-new_y][new_x] == '.':
                            # Make move
                            new_state = copy.deepcopy(current_state)
                            new_state[2-j][i] = '.'
                            new_state[2-new_y][new_x] = piece
                            
                            state_str = state_to_string(new_state)
                            if state_str not in seen_states:
                                seen_states.add(state_str)
                                move = f"{piece},{pos_to_notation((i,j))},{pos_to_notation((new_x,new_y))}"
                                queue.append((new_state, moves + [move], not is_black_turn))
    
    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")