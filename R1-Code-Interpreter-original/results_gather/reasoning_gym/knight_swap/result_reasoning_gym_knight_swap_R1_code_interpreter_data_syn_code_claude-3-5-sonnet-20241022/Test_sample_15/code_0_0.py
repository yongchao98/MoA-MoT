from collections import deque
import copy

def get_knight_moves(pos):
    row, col = pos
    moves = [
        (row+2, col+1), (row+2, col-1),
        (row-2, col+1), (row-2, col-1),
        (row+1, col+2), (row+1, col-2),
        (row-1, col+2), (row-1, col-2)
    ]
    return [(r,c) for r,c in moves if 0 <= r < 3 and 0 <= c < 4]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_final_state(board):
    # Check if white knights are in original black positions and vice versa
    white_positions = [(r,c) for r in range(3) for c in range(4) if board[r][c] == 'w']
    black_positions = [(r,c) for r in range(3) for c in range(4) if board[r][c] == 'B']
    return set(white_positions) == {(0,2), (0,3)} and set(black_positions) == {(1,1), (1,2)}

def solve_knight_swap():
    # Initial board state
    initial_board = [
        ['.', '.', '.', '.'],
        ['.', 'w', 'w', '.'],
        ['.', '.', 'B', 'B']
    ]
    
    # BFS queue: (board, moves_list, is_white_turn)
    queue = deque([(initial_board, [], True)])
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves, white_turn = queue.popleft()
        
        # Find current player's knights
        piece = 'w' if white_turn else 'B'
        knights = [(r,c) for r in range(3) for c in range(4) if current_board[r][c] == piece]
        
        # Try all possible moves for current knights
        for kr, kc in knights:
            for new_r, new_c in get_knight_moves((kr, kc)):
                if current_board[new_r][new_c] == '.':
                    # Make move
                    new_board = copy.deepcopy(current_board)
                    new_board[kr][kc] = '.'
                    new_board[new_r][new_c] = piece
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        move = f"{piece},{chr(ord('A')+kc)}{3-kr},{chr(ord('A')+new_c)}{3-new_r}"
                        new_moves = moves + [move]
                        
                        if is_final_state(new_board):
                            return new_moves
                        
                        queue.append((new_board, new_moves, not white_turn))
    
    return "No"

# Solve and print result
result = solve_knight_swap()
print(result if result == "No" else result)