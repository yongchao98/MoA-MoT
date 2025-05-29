from collections import deque

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 3]

def pos_to_coord(pos):
    return (ord(pos[0]) - ord('A'), int(pos[1]) - 1)

def coord_to_pos(coord):
    return chr(coord[0] + ord('A')) + str(coord[1] + 1)

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_final_state(white_positions, black_positions):
    # Check if white knights are at original black positions (B3,C3)
    # and black knights are at original white positions (A2,D2)
    return (sorted(white_positions) == [(1,0), (2,0)] and  # B3,C3
            sorted(black_positions) == [(0,1), (3,1)])     # A2,D2

def solve_knight_swap():
    # Initial positions
    initial_white = {(0,1), (3,1)}  # A2,D2
    initial_black = {(1,0), (2,0)}  # B3,C3
    
    # Initialize board
    initial_board = [['.'] * 4 for _ in range(3)]
    for wx, wy in initial_white:
        initial_board[wy][wx] = 'w'
    for bx, by in initial_black:
        initial_board[by][bx] = 'B'
    
    # (board, white_positions, black_positions, moves, is_white_turn)
    queue = deque([(initial_board, initial_white, initial_black, [], True)])
    visited = {(board_to_string(initial_board), True)}
    
    while queue:
        board, white_pos, black_pos, moves, is_white_turn = queue.popleft()
        current_positions = white_pos if is_white_turn else black_pos
        
        for pos in current_positions:
            for new_x, new_y in get_knight_moves(pos):
                if board[new_y][new_x] == '.':  # Empty square
                    # Create new state
                    new_board = [row[:] for row in board]
                    new_board[pos[1]][pos[0]] = '.'
                    new_board[new_y][new_x] = 'w' if is_white_turn else 'B'
                    
                    # Update positions
                    new_white_pos = white_pos.copy()
                    new_black_pos = black_pos.copy()
                    if is_white_turn:
                        new_white_pos.remove(pos)
                        new_white_pos.add((new_x, new_y))
                    else:
                        new_black_pos.remove(pos)
                        new_black_pos.add((new_x, new_y))
                    
                    state = (board_to_string(new_board), not is_white_turn)
                    if state not in visited:
                        visited.add(state)
                        from_pos = coord_to_pos(pos)
                        to_pos = coord_to_pos((new_x, new_y))
                        new_moves = moves + [f"{'w' if is_white_turn else 'B'},{from_pos},{to_pos}"]
                        
                        if is_final_state(new_white_pos, new_black_pos):
                            return new_moves
                        
                        queue.append((new_board, new_white_pos, new_black_pos, 
                                    new_moves, not is_white_turn))
    
    return "No"

solution = solve_knight_swap()
print(solution)