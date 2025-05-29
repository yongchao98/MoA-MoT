from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0 <= x < 3 and 0 <= y < 4]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_goal_state(board):
    # Check if white knights are at A1, B1 and black knights at B3, A2
    white_positions = set()
    black_positions = set()
    for i in range(3):
        for j in range(4):
            if board[i][j] == 'w':
                white_positions.add((i,j))
            elif board[i][j] == 'B':
                black_positions.add((i,j))
    
    return (0,0) in white_positions and (0,1) in white_positions and \
           (2,1) in black_positions and (1,0) in black_positions

def evaluate_state(board, target_white, target_black):
    # Heuristic function to estimate distance to goal
    score = 0
    for i in range(3):
        for j in range(4):
            if board[i][j] == 'w':
                score += min(abs(i-tx) + abs(j-ty) for tx,ty in target_white)
            elif board[i][j] == 'B':
                score += min(abs(i-tx) + abs(j-ty) for tx,ty in target_black)
    return score

def find_solution():
    initial_board = [
        ['.', 'w', '.', '.'],
        ['w', '.', '.', '.'],
        ['B', 'B', '.', '.']
    ]
    
    target_white = {(0,0), (0,1)}  # A1, B1
    target_black = {(2,1), (1,0)}  # B3, A2
    
    queue = deque([(initial_board, [], True, evaluate_state(initial_board, target_white, target_black))])
    visited = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn, _ = queue.popleft()
        
        if is_goal_state(board):
            return moves
            
        piece = 'w' if is_white_turn else 'B'
        pieces = [(i,j) for i in range(3) for j in range(4) if board[i][j] == piece]
        
        next_states = []
        for x, y in pieces:
            for new_x, new_y in get_knight_moves((x,y)):
                if board[new_x][new_y] == '.':
                    new_board = copy.deepcopy(board)
                    new_board[x][y] = '.'
                    new_board[new_x][new_y] = piece
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        move = f"{piece},{chr(y+65)}{3-x},{chr(new_y+65)}{3-new_x}"
                        new_moves = moves + [move]
                        score = evaluate_state(new_board, target_white, target_black)
                        next_states.append((new_board, new_moves, not is_white_turn, score))
        
        # Sort next states by score before adding to queue
        next_states.sort(key=lambda x: x[3])
        queue.extend(next_states)
    
    return None

solution = find_solution()
if solution:
    print(solution)
else:
    print("No")