from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [
        (x+2, y+1), (x+2, y-1), (x-2, y+1), (x-2, y-1),
        (x+1, y+2), (x+1, y-2), (x-1, y+2), (x-1, y-2)
    ]
    return [(x, y) for x, y in moves if 0 <= x < 4 and 0 <= y < 3]

def pos_to_notation(pos):
    return chr(pos[0] + ord('A')) + str(pos[1] + 1)

def state_to_string(state):
    return ''.join(''.join(row) for row in state)

def is_target_state(state):
    # Check if white knights are at original black positions and vice versa
    target_positions = {(0,0): 'w', (2,2): 'w', (1,0): 'B', (2,0): 'B'}
    for pos, piece in target_positions.items():
        if state[pos[1]][pos[0]] != piece:
            return False
    return True

def solve_knight_swap():
    # Initial state
    initial_state = [
        ['.', '.', 'B', '.'],
        ['.', '.', '.', '.'],
        ['B', 'w', 'w', '.']
    ]
    
    queue = deque([(initial_state, [], True)])  # state, moves, white_turn
    visited = {state_to_string(initial_state)}
    
    while queue:
        current_state, moves, white_turn = queue.popleft()
        
        if is_target_state(current_state):
            return moves
        
        # Find all pieces that can move
        piece = 'w' if white_turn else 'B'
        for y in range(3):
            for x in range(4):
                if current_state[y][x] == piece:
                    # Try all possible moves for this piece
                    for new_x, new_y in get_knight_moves((x, y)):
                        if current_state[new_y][new_x] == '.':
                            # Make move
                            new_state = copy.deepcopy(current_state)
                            new_state[y][x] = '.'
                            new_state[new_y][new_x] = piece
                            
                            state_str = state_to_string(new_state)
                            if state_str not in visited:
                                visited.add(state_str)
                                move = f"{piece},{pos_to_notation((x,y))},{pos_to_notation((new_x,new_y))}"
                                queue.append((new_state, moves + [move], not white_turn))
    
    return None

solution = solve_knight_swap()
if solution is None:
    print("<<<No>>>")
else:
    print(f"<<<{solution}>>>")