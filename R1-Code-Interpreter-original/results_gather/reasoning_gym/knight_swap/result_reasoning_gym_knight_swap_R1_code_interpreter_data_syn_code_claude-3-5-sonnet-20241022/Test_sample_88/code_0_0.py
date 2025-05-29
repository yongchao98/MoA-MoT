from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [
        (x+2, y+1), (x+2, y-1), (x-2, y+1), (x-2, y-1),
        (x+1, y+2), (x+1, y-2), (x-1, y+2), (x-1, y-2)
    ]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 3]

def pos_to_notation(pos):
    return chr(pos[0] + ord('A')) + str(pos[1] + 1)

def notation_to_pos(notation):
    return (ord(notation[0]) - ord('A'), int(notation[1]) - 1)

def is_target_state(board):
    # Check if white knights are at B1 and D1, black knights at A2 and C3
    return (board[0][1] == 'w' and board[0][3] == 'w' and 
            board[1][0] == 'B' and board[2][2] == 'B')

def solve_knight_swap():
    # Initial board state (0,0 is A1)
    initial_board = [
        ['.', 'B', '.', 'B'],  # Row 1
        ['w', '.', '.', '.'],  # Row 2
        ['.', '.', 'w', '.']   # Row 3
    ]
    
    visited = set()
    queue = deque([(initial_board, [], True)])  # board, moves, is_black_turn
    
    while queue:
        board, moves, is_black_turn = queue.popleft()
        board_tuple = tuple(tuple(row) for row in board)
        
        if board_tuple in visited:
            continue
        visited.add(board_tuple)
        
        if is_target_state(board):
            return moves
        
        # Find all knights of current color
        current_piece = 'B' if is_black_turn else 'w'
        knights = []
        for y in range(3):
            for x in range(4):
                if board[y][x] == current_piece:
                    knights.append((x, y))
        
        # Try all possible moves for each knight
        for knight in knights:
            for new_pos in get_knight_moves(knight):
                if board[new_pos[1]][new_pos[0]] == '.':
                    # Make move
                    new_board = copy.deepcopy(board)
                    new_board[knight[1]][knight[0]] = '.'
                    new_board[new_pos[1]][new_pos[0]] = current_piece
                    
                    # Record move
                    move = f"{current_piece},{pos_to_notation(knight)},{pos_to_notation(new_pos)}"
                    new_moves = moves + [move]
                    
                    queue.append((new_board, new_moves, not is_black_turn))
    
    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")