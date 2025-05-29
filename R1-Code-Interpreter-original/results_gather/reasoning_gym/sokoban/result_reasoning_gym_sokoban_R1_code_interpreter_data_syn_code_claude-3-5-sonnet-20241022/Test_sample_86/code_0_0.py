from collections import deque
import copy

# Initial state
initial_state = [
    "+++++++" ,
    "+-X@-+",
    "+--@--+",
    "X@$--+",
    "+----+",
    "+-*@XX+",
    "+----+",
    "+++++++"
]

class State:
    def __init__(self, board, player_pos, boxes, goals):
        self.board = board
        self.player_pos = player_pos
        self.boxes = boxes
        self.goals = goals
        
def get_initial_state(board):
    player_pos = None
    boxes = set()
    goals = set()
    
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] in ['*', '%']:  # player
                player_pos = (i, j)
            if board[i][j] in ['@', '$']:  # box
                boxes.add((i, j))
            if board[i][j] in ['X', '$', '%']:  # goal
                goals.add((i, j))
                
    return State(board, player_pos, boxes, goals)

def is_valid_move(state, dx, dy):
    new_x = state.player_pos[0] + dx
    new_y = state.player_pos[1] + dy
    
    # Check wall
    if state.board[new_x][new_y] == '+':
        return False
        
    # Check if pushing box
    if (new_x, new_y) in state.boxes:
        box_new_x = new_x + dx
        box_new_y = new_y + dy
        
        # Check if box can be pushed
        if state.board[box_new_x][box_new_y] == '+' or (box_new_x, box_new_y) in state.boxes:
            return False
            
    return True

def make_move(state, dx, dy, direction):
    new_state = State(state.board, state.player_pos, state.boxes.copy(), state.goals)
    
    new_x = state.player_pos[0] + dx
    new_y = state.player_pos[1] + dy
    
    # Move box if necessary
    if (new_x, new_y) in new_state.boxes:
        box_new_x = new_x + dx
        box_new_y = new_y + dy
        new_state.boxes.remove((new_x, new_y))
        new_state.boxes.add((box_new_x, box_new_y))
    
    new_state.player_pos = (new_x, new_y)
    return new_state, direction

def is_solved(state):
    return state.boxes == state.goals

def solve_sokoban(initial_board):
    moves = [(0, -1, 'L'), (0, 1, 'R'), (-1, 0, 'U'), (1, 0, 'D')]
    initial = get_initial_state(initial_board)
    
    queue = deque([(initial, "")])
    seen = set()
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_solved(current_state):
            return path
            
        state_hash = (current_state.player_pos, frozenset(current_state.boxes))
        if state_hash in seen:
            continue
        seen.add(state_hash)
        
        for dx, dy, direction in moves:
            if is_valid_move(current_state, dx, dy):
                new_state, move = make_move(current_state, dx, dy, direction)
                queue.append((new_state, path + move))
    
    return None

solution = solve_sokoban(initial_state)
print(solution)