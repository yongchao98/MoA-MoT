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

def pos_to_notation(pos):
    return f"{chr(65+pos[1])}{pos[0]+1}"

def solve_knight_swap():
    # Initial board state: 0=empty, 1=white knight, 2=black knight
    initial_board = [
        [0, 0, 0, 0],
        [0, 0, 1, 2],
        [2, 0, 0, 1]
    ]
    
    # Target board state (swapped positions)
    target_board = [
        [0, 0, 0, 0],
        [0, 0, 2, 1],
        [1, 0, 0, 2]
    ]
    
    # Find knight positions
    def get_knight_positions(board):
        white, black = [], []
        for i in range(3):
            for j in range(4):
                if board[i][j] == 1:
                    white.append((i,j))
                elif board[i][j] == 2:
                    black.append((i,j))
        return white, black
    
    # State representation: (board, white_to_move)
    initial_state = (tuple(map(tuple, initial_board)), True)
    visited = {initial_state}
    queue = deque([(initial_state, [])])
    
    while queue:
        current_state, moves = queue.popleft()
        current_board, white_turn = current_state
        current_board = [list(row) for row in current_board]
        
        # Check if we reached target state
        if current_board == target_board:
            return moves
        
        # Get current positions
        white_knights, black_knights = get_knight_positions(current_board)
        current_knights = white_knights if white_turn else black_knights
        
        # Try all possible moves for current color
        for knight in current_knights:
            for new_pos in get_knight_moves(knight):
                if current_board[new_pos[0]][new_pos[1]] == 0:  # If square is empty
                    # Make move
                    new_board = copy.deepcopy(current_board)
                    new_board[new_pos[0]][new_pos[1]] = new_board[knight[0]][knight[1]]
                    new_board[knight[0]][knight[1]] = 0
                    
                    new_state = (tuple(map(tuple, new_board)), not white_turn)
                    if new_state not in visited:
                        visited.add(new_state)
                        move = f"{'w' if white_turn else 'B'},{pos_to_notation(knight)},{pos_to_notation(new_pos)}"
                        queue.append((new_state, moves + [move]))
    
    return None

# Solve and print result
solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")