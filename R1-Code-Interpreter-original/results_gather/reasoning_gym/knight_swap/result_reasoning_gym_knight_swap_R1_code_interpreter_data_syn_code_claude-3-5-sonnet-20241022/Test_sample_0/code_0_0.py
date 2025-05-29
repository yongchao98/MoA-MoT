def get_valid_moves(board, pos):
    moves = [
        (-2, -1), (-2, 1), (-1, -2), (-1, 2),
        (1, -2), (1, 2), (2, -1), (2, 1)
    ]
    valid = []
    x, y = pos
    for dx, dy in moves:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < 4 and 0 <= new_y < 3 and board[new_y][new_x] == '.':
            valid.append((new_x, new_y))
    return valid

def pos_to_notation(pos):
    return chr(65 + pos[0]) + str(pos[1] + 1)

def state_to_board(state):
    board = [['.'] * 4 for _ in range(3)]
    for piece, pos in state.items():
        if pos:  # If piece hasn't reached target
            board[pos[1]][pos[0]] = piece[0]
    return board

def is_target_reached(state, targets):
    return all(state[piece] == target for piece, target in targets.items())

def solve_knight_swap():
    # Initial positions (x, y) - A1 is (0,0)
    initial_state = {
        'w1': (1,2),  # B3
        'w2': (2,0),  # C1
        'B1': (2,2),  # C3
        'B2': (0,1)   # A2
    }
    
    # Target positions
    targets = {
        'w1': (0,1),  # A2
        'w2': (2,2),  # C3
        'B1': (1,2),  # B3
        'B2': (2,0)   # C1
    }
    
    # BFS queue: (state, moves_list, white_turn)
    queue = [(initial_state, [], True)]
    seen = {str(initial_state)}
    
    while queue:
        state, moves, white_turn = queue.pop(0)
        
        if is_target_reached(state, targets):
            return moves
        
        # Get pieces that can move (white or black)
        pieces = [p for p in state.keys() if p.startswith('w' if white_turn else 'B')]
        board = state_to_board(state)
        
        for piece in pieces:
            if not state[piece]:  # Skip if piece has reached target
                continue
            
            curr_pos = state[piece]
            valid_moves = get_valid_moves(board, curr_pos)
            
            for new_pos in valid_moves:
                new_state = state.copy()
                new_state[piece] = new_pos
                
                state_key = str(new_state)
                if state_key not in seen:
                    seen.add(state_key)
                    move = f"{piece[0]},{pos_to_notation(curr_pos)},{pos_to_notation(new_pos)}"
                    queue.append((new_state, moves + [move], not white_turn))
    
    return None

# Solve and print result
result = solve_knight_swap()
if result:
    print(result)
else:
    print("No")