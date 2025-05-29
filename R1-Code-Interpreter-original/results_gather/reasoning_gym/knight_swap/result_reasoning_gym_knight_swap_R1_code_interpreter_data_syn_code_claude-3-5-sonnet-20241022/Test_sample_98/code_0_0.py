from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [
        (x+2, y+1), (x+2, y-1), (x-2, y+1), (x-2, y-1),
        (x+1, y+2), (x+1, y-2), (x-1, y+2), (x-1, y-2)
    ]
    return [(x, y) for x, y in moves if 0 <= x < 4 and 0 <= y < 3]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_goal_state(board):
    return (board[0][2] == 'B' and board[2][0] == 'w' and 
            board[2][1] == 'w' and board[2][2] == 'B')

def get_position_name(pos):
    x, y = pos
    return f"{chr(65+x)}{y+1}"

def solve_knight_swap():
    # Initial board state (0,0 is A1)
    initial_board = [
        ['.', '.', 'w', '.'],
        ['.', '.', '.', '.'],
        ['B', 'B', 'w', '.']
    ]
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    visited = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        if is_goal_state(board):
            return moves
        
        # Find current player's knights
        piece = 'w' if is_white_turn else 'B'
        knights = []
        for y in range(3):
            for x in range(4):
                if board[y][x] == piece:
                    knights.append((x, y))
        
        # Try all possible moves for each knight
        for knight in knights:
            for new_pos in get_knight_moves(knight):
                if board[new_pos[1]][new_pos[0]] == '.':
                    # Make move
                    new_board = copy.deepcopy(board)
                    new_board[knight[1]][knight[0]] = '.'
                    new_board[new_pos[1]][new_pos[0]] = piece
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        move = f"{piece},{get_position_name(knight)},{get_position_name(new_pos)}"
                        queue.append((new_board, moves + [move], not is_white_turn))
    
    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")