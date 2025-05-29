from collections import deque

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0 <= x < 3 and 0 <= y < 4]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_goal_state(board):
    # White knights should be at A1, B1; Black knights at B3, A2
    white_at_target = (board[2][0] == 'w' and board[2][1] == 'w')
    black_at_target = (board[0][1] == 'B' and board[1][0] == 'B')
    return white_at_target and black_at_target

def find_solution_with_depth(max_depth):
    initial_board = [
        ['.', 'w', '.', '.'],
        ['w', '.', '.', '.'],
        ['B', 'B', '.', '.']
    ]
    
    def dfs(board, moves, is_white_turn, depth, visited):
        if depth >= max_depth:
            return None
            
        if is_goal_state(board):
            return moves
            
        piece = 'w' if is_white_turn else 'B'
        
        # Find all pieces of current player
        pieces = []
        for i in range(3):
            for j in range(4):
                if board[i][j] == piece:
                    pieces.append((i,j))
        
        # Try all possible moves for each piece
        for x, y in pieces:
            for new_x, new_y in get_knight_moves((x,y)):
                if board[new_x][new_y] == '.':
                    # Make move
                    new_board = [row[:] for row in board]
                    new_board[x][y] = '.'
                    new_board[new_x][new_y] = piece
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        move = f"{piece},{chr(y+65)}{3-x},{chr(new_y+65)}{3-new_x}"
                        result = dfs(new_board, moves + [move], not is_white_turn, depth + 1, visited)
                        if result:
                            return result
                        visited.remove(board_str)
    
    return dfs(initial_board, [], True, 0, {board_to_string(initial_board)})

# Try increasing depths until solution is found
for depth in range(6, 11):  # We know solution exists within 10 moves
    solution = find_solution_with_depth(depth)
    if solution:
        print(solution)
        break