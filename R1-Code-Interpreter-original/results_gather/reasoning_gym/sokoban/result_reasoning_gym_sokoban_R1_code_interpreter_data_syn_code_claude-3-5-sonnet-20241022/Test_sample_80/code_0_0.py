from collections import deque
import heapq
import copy

class State:
    def __init__(self, board, player_pos, boxes, goals, parent=None, move=''):
        self.board = board
        self.player_pos = player_pos
        self.boxes = frozenset(boxes)
        self.goals = frozenset(goals)
        self.parent = parent
        self.move = move
        self.g = 0 if parent is None else parent.g + 1
        self.h = self.calculate_heuristic()
        
    def calculate_heuristic(self):
        total = 0
        for box in self.boxes:
            min_dist = float('inf')
            for goal in self.goals:
                dist = abs(box[0] - goal[0]) + abs(box[1] - goal[1])
                min_dist = min(min_dist, dist)
            total += min_dist
        return total
    
    def __lt__(self, other):
        return (self.g + self.h) < (other.g + other.h)
    
    def is_goal(self):
        return self.boxes == self.goals

def get_moves(state):
    moves = []
    directions = [('U', -1, 0), ('D', 1, 0), ('L', 0, -1), ('R', 0, 1)]
    
    for move, dy, dx in directions:
        new_y = state.player_pos[0] + dy
        new_x = state.player_pos[1] + dx
        
        if state.board[new_y][new_x] == '+':
            continue
            
        if (new_y, new_x) in state.boxes:
            box_new_y = new_y + dy
            box_new_x = new_x + dx
            
            if state.board[box_new_y][box_new_x] == '+' or (box_new_y, box_new_x) in state.boxes:
                continue
                
            new_boxes = set(state.boxes)
            new_boxes.remove((new_y, new_x))
            new_boxes.add((box_new_y, box_new_x))
            
            moves.append((move, (new_y, new_x), frozenset(new_boxes)))
        else:
            moves.append((move, (new_y, new_x), state.boxes))
    
    return moves

def solve_sokoban(board):
    # Parse initial state
    player_pos = None
    boxes = set()
    goals = set()
    
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] in ['*', '%']:
                player_pos = (i, j)
            if board[i][j] in ['@', '$']:
                boxes.add((i, j))
            if board[i][j] in ['X', '$', '%']:
                goals.add((i, j))
    
    initial_state = State(board, player_pos, boxes, goals)
    visited = set()
    pq = [(initial_state.g + initial_state.h, initial_state)]
    
    while pq:
        _, current = heapq.heappop(pq)
        
        if current.is_goal():
            # Reconstruct path
            path = []
            while current.parent:
                path.append(current.move)
                current = current.parent
            return ''.join(reversed(path))
        
        state_hash = (current.player_pos, current.boxes)
        if state_hash in visited:
            continue
        visited.add(state_hash)
        
        for move, new_player_pos, new_boxes in get_moves(current):
            new_state = State(board, new_player_pos, new_boxes, goals, current, move)
            heapq.heappush(pq, (new_state.g + new_state.h, new_state))
    
    return None

# Initialize the board
board = [
    list("+++++++++"),
    list("+-@X--+"),
    list("+-X-@--+"),
    list("+-@@*-X-+"),
    list("+X-$+--+"),
    list("+++++++++")
]

solution = solve_sokoban(board)
print(solution)