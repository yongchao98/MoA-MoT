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
    return [(r,c) for r,c in moves if 0 <= r < 3 and 0 <= c < 3]

def pos_to_notation(pos):
    return f"{chr(65+pos[1])}{3-pos[0]}"

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_target_reached(board):
    # Check if white knights are at A1, C1 and black knights are at C3, A2
    return (board[2][0] == 'w' and board[2][2] == 'w' and 
            board[0][2] == 'B' and board[1][0] == 'B')

def is_better_state(board):
    # Heuristic to guide search towards better states
    score = 0
    # Reward white knights being closer to bottom row
    for i in range(3):
        for j in range(3):
            if board[i][j] == 'w':
                score += i  # Higher row number (closer to bottom) is better
            elif board[i][j] == 'B':
                score += (2-i)  # Lower row number (closer to top) is better
    return score

def find_solution():
    initial_board = [
        ['.', '.', 'w'],  # row 3
        ['w', '.', '.'],  # row 2
        ['B', '.', 'B']   # row 1
    ]
    
    queue = deque([(initial_board, [], True)])
    seen = {board_to_string(initial_board)}
    best_solution = None
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        if is_target_reached(board):
            if best_solution is None or len(moves) < len(best_solution):
                best_solution = moves
                continue
        
        if best_solution and len(moves) >= len(best_solution):
            continue
            
        piece = 'w' if is_white_turn else 'B'
        pieces = []
        for i in range(3):
            for j in range(3):
                if board[i][j] == piece:
                    pieces.append((i,j))
        
        # Sort moves by heuristic score
        moves_with_scores = []
        for start_pos in pieces:
            for end_pos in get_knight_moves(start_pos):
                if board[end_pos[0]][end_pos[1]] == '.':
                    new_board = copy.deepcopy(board)
                    new_board[end_pos[0]][end_pos[1]] = piece
                    new_board[start_pos[0]][start_pos[1]] = '.'
                    score = is_better_state(new_board)
                    moves_with_scores.append((score, start_pos, end_pos))
        
        # Try moves in order of heuristic score
        for _, start_pos, end_pos in sorted(moves_with_scores, reverse=True):
            new_board = copy.deepcopy(board)
            new_board[end_pos[0]][end_pos[1]] = piece
            new_board[start_pos[0]][start_pos[1]] = '.'
            
            board_str = board_to_string(new_board)
            if board_str not in seen:
                seen.add(board_str)
                move = f"{piece},{pos_to_notation(start_pos)},{pos_to_notation(end_pos)}"
                queue.append((new_board, moves + [move], not is_white_turn))
    
    return best_solution

solution = find_solution()
if solution is None:
    print("<<<No>>>")
else:
    print(f"<<<{solution}>>>")