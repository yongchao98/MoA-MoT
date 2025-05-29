from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [
        (x+2, y+1), (x+2, y-1),
        (x-2, y+1), (x-2, y-1),
        (x+1, y+2), (x+1, y-2),
        (x-1, y+2), (x-1, y-2)
    ]
    return [(x, y) for x, y in moves if 0 <= x < 4 and 0 <= y < 3]

def pos_to_notation(pos):
    return chr(pos[0] + ord('A')) + str(pos[1] + 1)

def board_to_string(board):
    return ','.join(''.join(row) for row in board)

def is_target_reached(board):
    # Check if white knights are at B2 and D1, black knights at D3 and A2
    return (board[2][3] == 'B' and board[1][0] == 'B' and 
            board[2][1] == 'w' and board[0][3] == 'w')

def find_solution():
    initial_board = [
        ['.', '.', '.', 'w'],
        ['w', 'B', '.', '.'],
        ['.', '.', '.', 'B']
    ]
    
    # State: (board, is_white_turn, moves_list)
    start_state = (initial_board, True, [])
    visited = set()
    queue = deque([start_state])
    
    while queue:
        board, is_white_turn, moves = queue.popleft()
        board_str = board_to_string(board)
        if (board_str, is_white_turn) in visited:
            continue
        visited.add((board_str, is_white_turn))
        
        if is_target_reached(board):
            return moves
        
        # Find pieces to move
        piece = 'w' if is_white_turn else 'B'
        for i in range(3):
            for j in range(4):
                if board[i][j] == piece:
                    from_pos = (j, i)
                    # Try all possible moves
                    for to_x, to_y in get_knight_moves((j, i)):
                        if board[to_y][to_x] == '.':
                            new_board = copy.deepcopy(board)
                            new_board[i][j] = '.'
                            new_board[to_y][to_x] = piece
                            move = f"{piece},{pos_to_notation(from_pos)},{pos_to_notation((to_x, to_y))}"
                            new_moves = moves + [move]
                            queue.append((new_board, not is_white_turn, new_moves))
    
    return None

solution = find_solution()
if solution:
    print(solution)
else:
    print("No")