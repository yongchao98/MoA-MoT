from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [
        (x+2, y+1), (x+2, y-1), (x-2, y+1), (x-2, y-1),
        (x+1, y+2), (x+1, y-2), (x-1, y+2), (x-1, y-2)
    ]
    return [(x,y) for x,y in moves if 0 <= x < 3 and 0 <= y < 3]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def get_positions(board, piece):
    return [(i,j) for i in range(3) for j in range(3) if board[i][j] == piece]

def is_target_state(board):
    white_pos = get_positions(board, 'w')
    black_pos = get_positions(board, 'B')
    return sorted(white_pos) == [(1,0), (1,1)] and sorted(black_pos) == [(2,0), (0,0)]

def solve_knight_swap():
    initial_board = [
        ['w', '.', '.'],
        ['B', '.', '.'],
        ['w', 'B', '.']
    ]
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    visited = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        if is_target_state(board):
            return moves
            
        piece = 'w' if is_white_turn else 'B'
        positions = get_positions(board, piece)
        
        for pos in positions:
            for new_pos in get_knight_moves(pos):
                i, j = new_pos
                if board[i][j] == '.':
                    new_board = copy.deepcopy(board)
                    old_i, old_j = pos
                    new_board[i][j] = piece
                    new_board[old_i][old_j] = '.'
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        move = f"{piece},{chr(old_j+65)}{old_i+1},{chr(j+65)}{i+1}"
                        queue.append((new_board, moves + [move], not is_white_turn))
    
    return None

result = solve_knight_swap()
print(result if result else "No")