from collections import deque

def is_solved(boxes, goals):
    return all(box in goals for box in boxes)

def get_new_position(pos, move):
    if move == 'R':
        return (pos[0], pos[1] + 1)
    elif move == 'L':
        return (pos[0], pos[1] - 1)
    elif move == 'U':
        return (pos[0] - 1, pos[1])
    elif move == 'D':
        return (pos[0] + 1, pos[1])

def solve_sokoban_bfs():
    # Initial positions
    player_pos = (4, 5)
    boxes = {(2, 6), (3, 3), (4, 4), (5, 4), (6, 5)}
    goals = {(2, 3), (2, 5), (4, 7), (6, 2)}
    walls = {(0, i) for i in range(9)} | {(8, i) for i in range(9)} | {(i, 0) for i in range(9)} | {(i, 8) for i in range(9)}
    
    # BFS setup
    queue = deque([(player_pos, frozenset(boxes), "")])
    visited = set()
    
    while queue:
        player_pos, boxes, path = queue.popleft()
        
        if is_solved(boxes, goals):
            return path
        
        if (player_pos, boxes) in visited:
            continue
        visited.add((player_pos, boxes))
        
        for move in "RDLU":
            new_player_pos = get_new_position(player_pos, move)
            if new_player_pos in walls:
                continue
            
            new_boxes = set(boxes)
            if new_player_pos in boxes:
                new_box_pos = get_new_position(new_player_pos, move)
                if new_box_pos in walls or new_box_pos in boxes:
                    continue
                new_boxes.remove(new_player_pos)
                new_boxes.add(new_box_pos)
            
            queue.append((new_player_pos, frozenset(new_boxes), path + move))
    
    return None

solution = solve_sokoban_bfs()
print(solution)