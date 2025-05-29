from collections import deque

class State:
    def __init__(self, player, boxes, goals):
        self.player = player
        self.boxes = frozenset(boxes)  # immutable set for hashing
        self.goals = frozenset(goals)
        
    def is_complete(self):
        return self.boxes == self.goals
        
    def __hash__(self):
        return hash((self.player, self.boxes))
        
    def __eq__(self, other):
        return self.player == other.player and self.boxes == other.boxes

def parse_board(board):
    player = None
    boxes = set()
    goals = set()
    walls = set()
    
    for i in range(len(board)):
        for j in range(len(board[i])):
            pos = (i, j)
            cell = board[i][j]
            
            if cell == '+':
                walls.add(pos)
            elif cell == '@':
                player = pos
            elif cell == '$':
                boxes.add(pos)
            elif cell == 'X':
                goals.add(pos)
            elif cell == '*':
                player = pos
                goals.add(pos)
            elif cell == '%':
                player = pos
                goals.add(pos)
            
    return player, boxes, goals, walls

def get_valid_moves(state, walls):
    moves = []
    directions = [
        ('U', -1, 0),
        ('D', 1, 0),
        ('L', 0, -1),
        ('R', 0, 1)
    ]
    
    for direction, dx, dy in directions:
        new_player = (state.player[0] + dx, state.player[1] + dy)
        
        # Check if move is into wall
        if new_player in walls:
            continue
            
        # Moving to empty space
        if new_player not in state.boxes:
            new_state = State(new_player, state.boxes, state.goals)
            moves.append((new_state, direction))
            continue
            
        # Pushing box
        new_box = (new_player[0] + dx, new_player[1] + dy)
        if new_box not in walls and new_box not in state.boxes:
            new_boxes = set(state.boxes)
            new_boxes.remove(new_player)
            new_boxes.add(new_box)
            new_state = State(new_player, new_boxes, state.goals)
            moves.append((new_state, direction))
            
    return moves

def is_deadlock(state, walls):
    # Check for corner deadlocks
    for box in state.boxes:
        if box not in state.goals:
            x, y = box
            corners = [
                ((x+1, y), (x, y+1), (x+1, y+1)),
                ((x+1, y), (x, y-1), (x+1, y-1)),
                ((x-1, y), (x, y+1), (x-1, y+1)),
                ((x-1, y), (x, y-1), (x-1, y-1))
            ]
            
            for corner in corners:
                if all(pos in walls for pos in corner):
                    return True
    return False

def solve_sokoban(board):
    player, boxes, goals, walls = parse_board(board)
    initial_state = State(player, boxes, goals)
    
    if initial_state.is_complete():
        return ""
        
    queue = deque([(initial_state, "")])
    visited = {initial_state}
    
    while queue:
        current_state, path = queue.popleft()
        
        for next_state, move in get_valid_moves(current_state, walls):
            if next_state in visited:
                continue
                
            if is_deadlock(next_state, walls):
                continue
                
            if next_state.is_complete():
                return path + move
                
            visited.add(next_state)
            queue.append((next_state, path + move))
            
    return None

# Initial board
board = [
    ['+', '+', '+', '+', '+', '+', '+'],
    ['+', '+', '+', '+', '-', '$', '+'],
    ['+', '+', '-', '-', '-', '-', '+'],
    ['+', '-', '-', '-', '-', '-', '+'],
    ['+', 'X', '-', '-', '-', '-', '+'],
    ['+', '-', '-', '@', 'X', '-', '+'],
    ['+', '-', '-', '@', '+', '-', '+'],
    ['+', '$', '$', 'X', '-', '@', '+'],
    ['+', '+', '+', '+', '$', '*', '+'],
    ['+', '+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(board)
print(solution)