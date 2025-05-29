from collections import deque

# Define the initial state of the puzzle
initial_state = {
    'player': (5, 8),
    'boxes': {(2, 4), (2, 6), (4, 8), (5, 7)},
    'goals': {(1, 4), (3, 2), (3, 6), (4, 8)},
    'walls': {(0, i) for i in range(10)} | {(6, i) for i in range(10)} | {(i, 0) for i in range(7)} | {(i, 9) for i in range(7)} | {(2, 8), (2, 9), (3, 8), (3, 9), (4, 9), (5, 9)},
}

# Define possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1),
}

# Check if a position is free (not a wall or a box)
def is_free(position, boxes, walls):
    return position not in boxes and position not in walls

# Check if the puzzle is solved
def is_solved(boxes, goals):
    return boxes == goals

# Perform a breadth-first search to find the solution
def solve_sokoban(initial_state):
    queue = deque([(initial_state['player'], initial_state['boxes'], "")])
    visited = set()
    
    while queue:
        player, boxes, path = queue.popleft()
        
        if is_solved(boxes, initial_state['goals']):
            return path
        
        if (player, frozenset(boxes)) in visited:
            continue
        visited.add((player, frozenset(boxes)))
        
        for move, (dy, dx) in moves.items():
            new_player = (player[0] + dy, player[1] + dx)
            
            if is_free(new_player, boxes, initial_state['walls']):
                queue.append((new_player, boxes, path + move))
            elif new_player in boxes:
                new_box = (new_player[0] + dy, new_player[1] + dx)
                if is_free(new_box, boxes, initial_state['walls']):
                    new_boxes = set(boxes)
                    new_boxes.remove(new_player)
                    new_boxes.add(new_box)
                    queue.append((new_player, new_boxes, path + move))
    
    return None

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)