from heapq import heappop, heappush

def sokoban_solver():
    def manhattan_distance(p1, p2):
        return abs(p1[0] - p2[0]) + abs(p1[1] - p2[1])

    def heuristic(boxes, goals):
        return sum(min(manhattan_distance(box, goal) for goal in goals) for box in boxes)

    def is_deadlock(box, boxes, walls):
        # Check for simple deadlocks: box in a corner or against a wall
        x, y = box
        if ((x - 1, y) in walls or (x + 1, y) in walls) and ((x, y - 1) in walls or (x, y + 1) in walls):
            return True
        return False

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
    
    def search(path, g, threshold, player, boxes):
        f = g + heuristic(boxes, initial_state['goals'])
        if f > threshold:
            return f, None
        if boxes == initial_state['goals']:
            return f, path
        
        min_threshold = float('inf')
        for move, (dx, dy) in moves.items():
            new_player = (player[0] + dx, player[1] + dy)
            
            if new_player in initial_state['walls']:
                continue
            
            if new_player in boxes:
                new_box = (new_player[0] + dx, new_player[1] + dy)
                if new_box in initial_state['walls'] or new_box in boxes:
                    continue
                if is_deadlock(new_box, boxes, initial_state['walls']):
                    continue
                new_boxes = set(boxes)
                new_boxes.remove(new_player)
                new_boxes.add(new_box)
                new_boxes = frozenset(new_boxes)
            else:
                new_boxes = boxes
            
            if (new_player, new_boxes) not in visited:
                visited.add((new_player, new_boxes))
                t, result = search(path + move, g + 1, threshold, new_player, new_boxes)
                if result is not None:
                    return t, result
                if t < min_threshold:
                    min_threshold = t
        return min_threshold, None
    
    threshold = heuristic(initial_state['boxes'], initial_state['goals'])
    visited = set()
    path = ""
    
    while True:
        t, result = search(path, 0, threshold, initial_state['player'], frozenset(initial_state['boxes']))
        if result is not None:
            return result
        if t == float('inf'):
            return "No solution found"
        threshold = t

# Execute the solver
solution = sokoban_solver()
print(solution)