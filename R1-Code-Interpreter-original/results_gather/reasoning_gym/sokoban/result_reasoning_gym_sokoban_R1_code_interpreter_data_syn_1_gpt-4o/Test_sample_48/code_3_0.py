from collections import deque

def is_valid_move(player_pos, box_pos, move, walls):
    # Calculate new player position
    if move == 'R':
        new_player_pos = (player_pos[0], player_pos[1] + 1)
    elif move == 'L':
        new_player_pos = (player_pos[0], player_pos[1] - 1)
    elif move == 'U':
        new_player_pos = (player_pos[0] - 1, player_pos[1])
    elif move == 'D':
        new_player_pos = (player_pos[0] + 1, player_pos[1])
    
    # Check if new player position is a wall
    if new_player_pos in walls:
        return False
    
    # Check if new player position is a box
    if new_player_pos in box_pos:
        # Calculate new box position
        if move == 'R':
            new_box_pos = (new_player_pos[0], new_player_pos[1] + 1)
        elif move == 'L':
            new_box_pos = (new_player_pos[0], new_player_pos[1] - 1)
        elif move == 'U':
            new_box_pos = (new_player_pos[0] - 1, new_player_pos[1])
        elif move == 'D':
            new_box_pos = (new_player_pos[0] + 1, new_player_pos[1])
        
        # Check if new box position is a wall or another box
        if new_box_pos in walls or new_box_pos in box_pos:
            return False
    
    return True

def solve_sokoban():
    # Initial positions
    player_pos = (7, 2)
    boxes = [(4, 2), (6, 1), (7, 3)]
    goals = [(2, 2), (3, 1), (7, 4)]
    walls = {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7),
             (1, 0), (1, 1), (1, 2), (1, 3), (1, 5), (1, 6), (1, 7),
             (2, 0), (2, 1), (2, 3), (2, 4), (2, 6), (2, 7),
             (3, 0), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7),
             (4, 0), (4, 3), (4, 4), (4, 5), (4, 6), (4, 7),
             (5, 0), (5, 5), (5, 6), (5, 7),
             (6, 0), (6, 6), (6, 7),
             (7, 0), (7, 5), (7, 6), (7, 7),
             (8, 0), (8, 1), (8, 2), (8, 3), (8, 4), (8, 5), (8, 6), (8, 7)}
    
    # BFS setup
    queue = deque([(player_pos, tuple(boxes), "")])
    visited = set()
    visited.add((player_pos, tuple(boxes)))
    
    # Possible moves
    directions = ['R', 'L', 'U', 'D']
    
    while queue:
        current_player_pos, current_boxes, path = queue.popleft()
        
        # Check if all boxes are on goals
        if all(box in goals for box in current_boxes):
            print("Solution found:", path)
            return
        
        for move in directions:
            if is_valid_move(current_player_pos, current_boxes, move, walls):
                # Calculate new player position
                if move == 'R':
                    new_player_pos = (current_player_pos[0], current_player_pos[1] + 1)
                elif move == 'L':
                    new_player_pos = (current_player_pos[0], current_player_pos[1] - 1)
                elif move == 'U':
                    new_player_pos = (current_player_pos[0] - 1, current_player_pos[1])
                elif move == 'D':
                    new_player_pos = (current_player_pos[0] + 1, current_player_pos[1])
                
                # Calculate new box positions
                new_boxes = list(current_boxes)
                if new_player_pos in current_boxes:
                    box_index = current_boxes.index(new_player_pos)
                    if move == 'R':
                        new_boxes[box_index] = (new_boxes[box_index][0], new_boxes[box_index][1] + 1)
                    elif move == 'L':
                        new_boxes[box_index] = (new_boxes[box_index][0], new_boxes[box_index][1] - 1)
                    elif move == 'U':
                        new_boxes[box_index] = (new_boxes[box_index][0] - 1, new_boxes[box_index][1])
                    elif move == 'D':
                        new_boxes[box_index] = (new_boxes[box_index][0] + 1, new_boxes[box_index][1])
                
                new_state = (new_player_pos, tuple(new_boxes))
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_player_pos, tuple(new_boxes), path + move))
    
    print("Solution not found")

solve_sokoban()