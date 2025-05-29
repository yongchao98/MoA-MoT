from collections import deque

def sokoban_solver():
    # Define the initial state
    initial_state = {
        'player': (3, 2),
        'boxes': {(2, 2), (2, 5), (2, 7), (5, 2), (5, 7), (6, 4), (6, 7)},
        'goals': {(3, 5), (4, 2), (5, 3), (6, 3), (6, 6)},
        'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8), (0, 9),
                  (1, 0), (1, 9), (2, 0), (2, 9), (3, 0), (3, 9), (4, 0), (4, 9), (5, 0), (5, 9),
                  (6, 0), (6, 9), (7, 0), (7, 1), (7, 2), (7, 3), (7, 9), (8, 0), (8, 1), (8, 2),
                  (8, 3), (8, 4), (8, 5), (8, 6), (8, 7), (8, 8), (8, 9)},
    }
    
    # Define possible moves
    moves = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    
    # BFS setup
    queue = deque([(initial_state['player'], frozenset(initial_state['boxes']), "")])
    visited = set()
    
    while queue:
        player, boxes, path = queue.popleft()
        
        # Check if all boxes are on goals
        if boxes == initial_state['goals']:
            return path
        
        # Explore possible moves
        for move, (dx, dy) in moves.items():
            new_player = (player[0] + dx, player[1] + dy)
            
            if new_player in initial_state['walls']:
                continue
            
            if new_player in boxes:
                new_box = (new_player[0] + dx, new_player[1] + dy)
                if new_box in initial_state['walls'] or new_box in boxes:
                    continue
                new_boxes = set(boxes)
                new_boxes.remove(new_player)
                new_boxes.add(new_box)
                new_boxes = frozenset(new_boxes)
            else:
                new_boxes = boxes
            
            if (new_player, new_boxes) not in visited:
                visited.add((new_player, new_boxes))
                queue.append((new_player, new_boxes, path + move))
    
    return "No solution found"

# Execute the solver
solution = sokoban_solver()
print(solution)