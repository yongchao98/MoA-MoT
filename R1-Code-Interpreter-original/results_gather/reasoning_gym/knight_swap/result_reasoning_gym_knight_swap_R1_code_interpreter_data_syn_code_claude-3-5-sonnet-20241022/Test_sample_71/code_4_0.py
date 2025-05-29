from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [
        (x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
        (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)
    ]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 3]

def pos_to_notation(pos):
    return chr(65 + pos[0]) + str(3 - pos[1])

def board_to_string(board):
    return tuple(tuple(row) for row in board)

def is_goal_state(board):
    return (board[0][0] == 'w' and board[0][3] == 'w' and 
            board[1][2] == 'B' and board[1][3] == 'B')

def get_piece_positions(board, piece):
    return [(x,y) for y in range(3) for x in range(4) if board[y][x] == piece]

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def heuristic(board):
    white_pos = get_piece_positions(board, 'w')
    black_pos = get_piece_positions(board, 'B')
    target_white = [(0,0), (0,3)]  # A3, D3
    target_black = [(2,2), (2,3)]  # C2, D2
    
    white_dist = sum(min(manhattan_distance(w, t) for t in target_white) for w in white_pos)
    black_dist = sum(min(manhattan_distance(b, t) for t in target_black) for b in black_pos)
    return white_dist + black_dist

def solve_puzzle():
    initial_board = [
        ['B', '.', '.', 'B'],
        ['.', '.', 'w', 'w'],
        ['.', '.', '.', '.']
    ]
    
    queue = [(0, initial_board, [], True)]  # priority, board, moves, is_white_turn
    seen = {(board_to_string(initial_board), True)}
    
    while queue:
        queue.sort()  # Sort by priority
        _, board, moves, is_white_turn = queue.pop(0)
        
        if is_goal_state(board):
            return moves
        
        piece = 'w' if is_white_turn else 'B'
        pieces = get_piece_positions(board, piece)
        
        for from_pos in pieces:
            for to_pos in get_knight_moves(from_pos):
                x2, y2 = to_pos
                if board[y2][x2] == '.':
                    new_board = copy.deepcopy(board)
                    new_board[from_pos[1]][from_pos[0]] = '.'
                    new_board[y2][x2] = piece
                    
                    state = (board_to_string(new_board), not is_white_turn)
                    if state not in seen:
                        seen.add(state)
                        move = f"{piece},{pos_to_notation(from_pos)},{pos_to_notation(to_pos)}"
                        new_moves = moves + [move]
                        priority = len(new_moves) + heuristic(new_board)
                        queue.append((priority, new_board, new_moves, not is_white_turn))
    
    return "No"

# Solve and print result
result = solve_puzzle()
print(result)