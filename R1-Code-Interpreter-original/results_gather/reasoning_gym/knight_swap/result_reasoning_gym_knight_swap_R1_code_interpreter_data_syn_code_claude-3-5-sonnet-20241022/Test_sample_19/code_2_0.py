from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [
        (x+2, y+1), (x+2, y-1), (x-2, y+1), (x-2, y-1),
        (x+1, y+2), (x+1, y-2), (x-1, y+2), (x-1, y-2)
    ]
    return [(x,y) for x,y in moves if 0 <= x < 4 and y >= 0 and y < 3]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def heuristic_score(board):
    # Find positions of all pieces
    white_pos = []
    black_pos = []
    target_white = [(2,2), (3,2)]  # C1, D1 positions
    target_black = [(1,0), (3,1)]  # B3, D2 positions
    
    for i in range(3):
        for j in range(4):
            if board[i][j] == 'w':
                white_pos.append((j,i))
            elif board[i][j] == 'B':
                black_pos.append((j,i))
    
    # Calculate total distance to targets
    score = 0
    for w in white_pos:
        score += min(manhattan_distance(w, t) for t in target_white)
    for b in black_pos:
        score += min(manhattan_distance(b, t) for t in target_black)
    return score

def is_target_state(board):
    return (board[2][2] == 'w' and board[2][3] == 'w' and 
            board[0][1] == 'B' and board[1][3] == 'B')

def find_solution():
    initial_board = [
        ['.', 'w', '.', '.'],  # row 3
        ['.', '.', '.', 'w'],  # row 2
        ['.', '.', 'B', 'B']   # row 1
    ]
    
    # Priority queue with (heuristic_score, board, moves, turn)
    queue = [(heuristic_score(initial_board), initial_board, [], True)]
    seen = {board_to_string(initial_board)}
    max_depth = 8  # Limit search depth
    
    while queue:
        _, board, moves, white_turn = queue.pop(0)
        
        if len(moves) >= max_depth:
            continue
            
        pieces = []
        for i in range(3):
            for j in range(4):
                if board[i][j] == ('w' if white_turn else 'B'):
                    pieces.append((j, i))
        
        next_states = []
        for piece in pieces:
            for new_pos in get_knight_moves(piece):
                if board[new_pos[1]][new_pos[0]] == '.':
                    new_board = copy.deepcopy(board)
                    new_board[piece[1]][piece[0]] = '.'
                    new_board[new_pos[1]][new_pos[0]] = 'w' if white_turn else 'B'
                    
                    board_str = board_to_string(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        
                        from_sq = f"{chr(65+piece[0])}{3-piece[1]}"
                        to_sq = f"{chr(65+new_pos[0])}{3-new_pos[1]}"
                        move = f"{'w' if white_turn else 'B'},{from_sq},{to_sq}"
                        new_moves = moves + [move]
                        
                        if is_target_state(new_board):
                            return new_moves
                        
                        score = heuristic_score(new_board)
                        next_states.append((score, new_board, new_moves, not white_turn))
        
        # Sort next states by heuristic score
        next_states.sort(key=lambda x: x[0])
        queue.extend(next_states)
    
    return "No"

result = find_solution()
print(result)