from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [
        (x+2, y+1), (x+2, y-1), (x-2, y+1), (x-2, y-1),
        (x+1, y+2), (x+1, y-2), (x-1, y+2), (x-1, y-2)
    ]
    return [(x, y) for x, y in moves if 0 <= x < 4 and 0 <= y < 3]

def pos_to_notation(pos):
    return chr(65 + pos[0]) + str(pos[1] + 1)

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_goal_state(board):
    # Check if white knights are at B3 and C2, black knights at A3 and D2
    return (board[2][1] == 'w' and board[1][2] == 'w' and 
            board[2][0] == 'B' and board[1][3] == 'B')

def solve_knight_swap():
    # Initial board state
    initial_board = [
        ['.','.','.','.'],
        ['.','.','B','w'],
        ['w','B','.','.']
    ]
    
    visited = set()
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        board_str = board_to_string(board)
        
        if board_str in visited:
            continue
        visited.add(board_str)
        
        if is_goal_state(board):
            return moves
        
        # Find all pieces of current color
        piece = 'w' if is_white_turn else 'B'
        pieces = []
        for y in range(3):
            for x in range(4):
                if board[y][x] == piece:
                    pieces.append((x, y))
        
        # Try all possible moves for each piece
        for from_pos in pieces:
            for to_pos in get_knight_moves(from_pos):
                x, y = to_pos
                if board[y][x] == '.':  # If target square is empty
                    new_board = copy.deepcopy(board)
                    # Make the move
                    new_board[from_pos[1]][from_pos[0]] = '.'
                    new_board[y][x] = piece
                    # Record the move
                    move = f"{piece},{pos_to_notation(from_pos)},{pos_to_notation(to_pos)}"
                    queue.append((new_board, moves + [move], not is_white_turn))
    
    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")