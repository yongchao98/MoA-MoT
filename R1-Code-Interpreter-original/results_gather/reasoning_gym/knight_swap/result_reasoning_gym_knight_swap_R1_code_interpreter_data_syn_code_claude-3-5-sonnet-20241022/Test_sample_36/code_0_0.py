from collections import deque

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 4]

def pos_to_coord(pos):
    return (ord(pos[0]) - ord('A'), int(pos[1]) - 1)

def coord_to_pos(coord):
    return chr(coord[0] + ord('A')) + str(coord[1] + 1)

def get_board_state(board):
    return ''.join(''.join(row) for row in board)

def is_target_state(board):
    return (board[2][0] == 'B' and board[2][1] == 'B' and
            board[1][0] == 'w' and board[1][3] == 'w')

def solve_knight_swap():
    # Initial board
    initial_board = [
        ['B', 'B', '.', '.'],
        ['w', '.', '.', 'w'],
        ['.', '.', '.', '.']
    ]
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    visited = {get_board_state(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        # Find all knights of current color
        knights = []
        for i in range(3):
            for j in range(4):
                if board[i][j] == ('w' if is_white_turn else 'B'):
                    knights.append((j, 2-i))  # Convert to coordinate system
        
        # Try all possible moves for each knight
        for knight in knights:
            for new_pos in get_knight_moves(knight):
                if board[2-new_pos[1]][new_pos[0]] == '.':  # Empty square
                    # Make move
                    new_board = [row[:] for row in board]
                    new_board[2-knight[1]][knight[0]] = '.'
                    new_board[2-new_pos[1]][new_pos[0]] = 'w' if is_white_turn else 'B'
                    
                    board_state = get_board_state(new_board)
                    if board_state not in visited:
                        visited.add(board_state)
                        from_pos = coord_to_pos((knight[0], 2-knight[1]))
                        to_pos = coord_to_pos((new_pos[0], 2-new_pos[1]))
                        new_moves = moves + [f"{'w' if is_white_turn else 'B'},{from_pos},{to_pos}"]
                        
                        if is_target_state(new_board):
                            return new_moves
                        
                        queue.append((new_board, new_moves, not is_white_turn))
    
    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")