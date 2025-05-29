from collections import deque

def get_next_pos(x, y, direction):
    if direction == 'U': return (x-1, y)
    if direction == 'D': return (x+1, y)
    if direction == 'L': return (x, y-1)
    if direction == 'R': return (x, y+1)

def solve_sokoban(board):
    height, width = len(board), len(board[0])
    
    # Find initial state
    player_pos = None
    boxes = set()
    goals = set()
    walls = set()
    
    for i in range(height):
        for j in range(width):
            if board[i][j] in ['@', '*']:
                player_pos = (i, j)
            if board[i][j] in ['$']:
                boxes.add((i, j))
            if board[i][j] in ['X', '*']:
                goals.add((i, j))
            if board[i][j] == '+':
                walls.add((i, j))
    
    def is_valid_move(pos):
        return pos not in walls and 0 <= pos[0] < height and 0 <= pos[1] < width
    
    def is_deadlock(box_pos, boxes):
        # Quick deadlock check for corners
        if (box_pos[0]-1, box_pos[1]) in walls and (box_pos[0], box_pos[1]-1) in walls:
            return True
        if (box_pos[0]-1, box_pos[1]) in walls and (box_pos[0], box_pos[1]+1) in walls:
            return True
        if (box_pos[0]+1, box_pos[1]) in walls and (box_pos[0], box_pos[1]-1) in walls:
            return True
        if (box_pos[0]+1, box_pos[1]) in walls and (box_pos[0], box_pos[1]+1) in walls:
            return True
        return False
    
    # BFS with state pruning
    queue = deque([(player_pos, frozenset(boxes), "")])
    visited = {(player_pos, frozenset(boxes))}
    max_moves = 30  # Limit search depth
    
    while queue:
        curr_player, curr_boxes, path = queue.popleft()
        
        if len(path) > max_moves:
            continue
            
        if curr_boxes == goals:
            return path
        
        for direction in ['U', 'D', 'L', 'R']:
            next_player = get_next_pos(*curr_player, direction)
            
            if not is_valid_move(next_player):
                continue
                
            new_boxes = set(curr_boxes)
            
            if next_player in curr_boxes:
                next_box = get_next_pos(*next_player, direction)
                
                if not is_valid_move(next_box) or next_box in curr_boxes:
                    continue
                
                if is_deadlock(next_box, new_boxes):
                    continue
                    
                new_boxes.remove(next_player)
                new_boxes.add(next_box)
            
            new_boxes = frozenset(new_boxes)
            new_state = (next_player, new_boxes)
            
            if new_state not in visited:
                visited.add(new_state)
                queue.append((next_player, new_boxes, path + direction))
    
    return None

# Initialize board
board = [
    ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+'],
    ['+', '-', '-', '-', '-', '-', '-', '-', '-', '+'],
    ['+', '-', 'X', '-', '$', '$', '-', '@', '-', '+'],
    ['+', '-', '-', '-', 'X', 'X', 'X', '-', '-', '+'],
    ['+', '+', 'X', '-', '@', '@', '@', '@', '*', '+'],
    ['+', '+', '-', '-', '-', '-', '-', '-', '-', '+'],
    ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+']
]

print(solve_sokoban(board))