def solve_sokoban():
    # Initial positions
    player_pos = (3, 4)
    boxes = [(2, 4), (5, 4), (6, 2), (6, 3)]
    goals = [(2, 2), (4, 3), (4, 4)]
    
    # New sequence of moves
    moves = "UULLLDRRRDDLUURRDDLL"
    
    # Simulate the moves
    for move in moves:
        if move == 'U':
            player_pos = (player_pos[0] - 1, player_pos[1])
        elif move == 'D':
            player_pos = (player_pos[0] + 1, player_pos[1])
        elif move == 'L':
            player_pos = (player_pos[0], player_pos[1] - 1)
        elif move == 'R':
            player_pos = (player_pos[0], player_pos[1] + 1)
        
        # Check if player is pushing a box
        if player_pos in boxes:
            box_index = boxes.index(player_pos)
            if move == 'U':
                boxes[box_index] = (boxes[box_index][0] - 1, boxes[box_index][1])
            elif move == 'D':
                boxes[box_index] = (boxes[box_index][0] + 1, boxes[box_index][1])
            elif move == 'L':
                boxes[box_index] = (boxes[box_index][0], boxes[box_index][1] - 1)
            elif move == 'R':
                boxes[box_index] = (boxes[box_index][0], boxes[box_index][1] + 1)
    
    # Check if all boxes are on goals
    boxes_on_goals = all(box in goals for box in boxes)
    return boxes_on_goals

print(solve_sokoban())