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
    return chr(65 + pos[0]) + str(pos[1] + 1)

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_goal_state(board):
    # Check if white knights are at black's initial positions and vice versa
    return (board[2][0] == 'w' and board[2][3] == 'w' and 
            board[1][2] == 'B' and board[1][3] == 'B')

def solve_puzzle():
    # Initial board state
    initial_board = [
        ['B', '.', '.', 'B'],
        ['.', '.', 'w', 'w'],
        ['.', '.', '.', '.']
    ]
    
    # Queue stores (board, moves_list, is_white_turn)
    queue = deque([(initial_board, [], True)])
    seen = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        piece = 'w' if is_white_turn else 'B'
        
        # Find all pieces of current player
        pieces = []
        for y in range(3):
            for x in range(4):
                if board[y][x] == piece:
                    pieces.append((x,y))
        
        # Try all possible moves for each piece
        for from_pos in pieces:
            for to_pos in get_knight_moves(from_pos):
                x2, y2 = to_pos
                if board[y2][x2] == '.':  # If target square is empty
                    # Make move
                    new_board = copy.deepcopy(board)
                    new_board[from_pos[1]][from_pos[0]] = '.'
                    new_board[y2][x2] = piece
                    
                    board_str = board_to_string(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        move = f"{piece},{pos_to_notation(from_pos)},{pos_to_notation(to_pos)}"
                        new_moves = moves + [move]
                        
                        if is_goal_state(new_board):
                            return new_moves
                        
                        queue.append((new_board, new_moves, not is_white_turn))
    
    return "No"

# Solve and print result
result = solve_puzzle()
print(result)