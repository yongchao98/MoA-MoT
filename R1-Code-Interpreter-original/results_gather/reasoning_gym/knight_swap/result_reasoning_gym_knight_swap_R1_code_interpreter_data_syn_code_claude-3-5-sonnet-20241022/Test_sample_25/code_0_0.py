from collections import deque

def get_knight_moves(pos):
    x, y = pos
    moves = [
        (x+2, y+1), (x+2, y-1), (x-2, y+1), (x-2, y-1),
        (x+1, y+2), (x+1, y-2), (x-1, y+2), (x-1, y-2)
    ]
    return [(x, y) for x, y in moves if 0 <= x < 4 and 0 <= y < 3]

def pos_to_notation(pos):
    return chr(65 + pos[0]) + str(pos[1] + 1)

def solve_knight_swap():
    # Initial positions (x, y) where A1 is (0,0)
    white_knights = {(1,2), (0,1)}  # B3, A2
    black_knights = {(3,1), (0,0)}  # D2, A1
    target_whites = {(0,0), (3,1)}  # A1, D2
    target_blacks = {(1,2), (0,1)}  # B3, A2
    
    initial_state = (frozenset(white_knights), frozenset(black_knights), True)  # True = white's turn
    visited = {initial_state}
    queue = deque([(initial_state, [])])
    
    while queue:
        state, moves = queue.popleft()
        whites, blacks, white_turn = state
        
        # Check if we reached target positions
        if whites == target_whites and blacks == target_blacks:
            return moves
        
        current_pieces = whites if white_turn else blacks
        other_pieces = blacks if white_turn else whites
        all_pieces = whites | blacks
        
        # Try all possible moves for current player
        for piece in current_pieces:
            for new_pos in get_knight_moves(piece):
                if new_pos not in all_pieces:  # Square is empty
                    new_current = frozenset(p if p != piece else new_pos for p in current_pieces)
                    new_state = (new_current, other_pieces) if white_turn else (other_pieces, new_current)
                    new_state = (*new_state, not white_turn)
                    
                    if new_state not in visited:
                        visited.add(new_state)
                        new_moves = moves + [f"{'w' if white_turn else 'B'},{pos_to_notation(piece)},{pos_to_notation(new_pos)}"]
                        queue.append((new_state, new_moves))
    
    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")