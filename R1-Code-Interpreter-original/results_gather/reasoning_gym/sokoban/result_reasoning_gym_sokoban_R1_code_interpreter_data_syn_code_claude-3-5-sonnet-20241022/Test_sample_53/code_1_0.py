from collections import deque

def get_player_and_boxes(board):
    player = None
    boxes = []
    goals = []
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] == '*':  # player
                player = (i, j)
            elif board[i][j] == '%':  # player on goal
                player = (i, j)
                goals.append((i, j))
            elif board[i][j] == '@':  # box
                boxes.append((i, j))
            elif board[i][j] == '$':  # box on goal
                boxes.append((i, j))
                goals.append((i, j))
            elif board[i][j] == 'X':  # goal
                goals.append((i, j))
    return player, set(boxes), set(goals)

def is_valid_move(pos, board):
    return (0 <= pos[0] < len(board) and 
            0 <= pos[1] < len(board[0]) and 
            board[pos[0]][pos[1]] != '+')

def try_move(player, boxes, direction, board):
    new_player = (player[0] + direction[0], player[1] + direction[1])
    
    if not is_valid_move(new_player, board):
        return None, None
        
    new_boxes = boxes.copy()
    
    # If moving into a box
    if new_player in boxes:
        new_box = (new_player[0] + direction[0], new_player[1] + direction[1])
        if not is_valid_move(new_box, board) or new_box in boxes:
            return None, None
        new_boxes.remove(new_player)
        new_boxes.add(new_box)
    
    return new_player, new_boxes

def solve_sokoban(board):
    player, boxes, goals = get_player_and_boxes(board)
    
    directions = {
        'U': (-1, 0),
        'D': (1, 0),
        'L': (0, -1),
        'R': (0, 1)
    }
    
    initial_state = (player, frozenset(boxes))
    queue = deque([(initial_state, "")])
    seen = {initial_state}
    
    while queue:
        (current_player, current_boxes), path = queue.popleft()
        
        # Check if solved
        if all(box in goals for box in current_boxes):
            return path
            
        # Try each direction
        for move, dir in directions.items():
            new_player, new_boxes = try_move(current_player, current_boxes, dir, board)
            
            if new_player is not None:
                new_state = (new_player, frozenset(new_boxes))
                if new_state not in seen:
                    seen.add(new_state)
                    queue.append((new_state, path + move))
                    
                    # Early success check
                    if all(box in goals for box in new_boxes):
                        return path + move
    
    return None

# Initialize the board
board = [
    list("+++++++" ),
    list("+++++X+" ),
    list("++++-@+" ),
    list("+++-@*+" ),
    list("++$-X$+" ),
    list("+-+-+++" ),
    list("+@+++++" ),
    list("+X+++++" ),
    list("+++$+$+" ),
    list("+++++++" )
]

# Print initial state for verification
print("Initial state:")
for row in board:
    print(''.join(row))
print("\nSolving...")

solution = solve_sokoban(board)
print("\nSolution:", solution if solution else "No solution found")