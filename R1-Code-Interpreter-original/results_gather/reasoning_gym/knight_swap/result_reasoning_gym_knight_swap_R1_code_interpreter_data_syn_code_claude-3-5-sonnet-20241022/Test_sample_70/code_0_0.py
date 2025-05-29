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

def is_target_state(board):
    # Check if white knights are where black knights started and vice versa
    return (board[2][3] == 'w' and board[0][3] == 'w' and 
            board[1][0] == 'B' and board[0][0] == 'B')

def find_solution():
    # Initial board state
    initial_board = [
        ['.', '.', '.', 'B'],
        ['w', '.', '.', '.'],
        ['w', '.', '.', 'B']
    ]
    
    # Queue for BFS: (board, moves_list, is_white_turn)
    queue = deque([(initial_board, [], True)])
    visited = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        # Find current knight positions
        knights = []
        for i in range(3):
            for j in range(4):
                if board[i][j] in ('w' if is_white_turn else 'B'):
                    knights.append((i, j))
        
        # Try all possible moves for current player's knights
        for knight_pos in knights:
            r, c = knight_pos
            piece = board[r][c]
            
            # Get all valid knight moves
            for new_r, new_c in get_knight_moves(knight_pos):
                if board[new_r][new_c] == '.':  # If target square is empty
                    # Make move
                    new_board = copy.deepcopy(board)
                    new_board[r][c] = '.'
                    new_board[new_r][new_c] = piece
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        move = f"{piece},{chr(65+c)}{r+1},{chr(65+new_c)}{new_r+1}"
                        new_moves = moves + [move]
                        
                        if is_target_state(new_board):
                            return new_moves
                        
                        queue.append((new_board, new_moves, not is_white_turn))
    
    return None

# Solve and print result
solution = find_solution()
if solution:
    print(solution)
else:
    print("No")