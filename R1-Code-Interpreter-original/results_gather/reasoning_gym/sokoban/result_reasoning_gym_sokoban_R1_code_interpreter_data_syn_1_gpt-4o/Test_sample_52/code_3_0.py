def solve_sokoban():
    # Initial positions
    player_pos = (4, 5)
    boxes = [(2, 3), (3, 4), (4, 4)]
    goals = [(3, 2), (4, 4)]
    
    # New sequence of moves
    moves = "ULLDRRUL"
    
    # Apply moves
    for move in moves:
        if move == 'L':
            player_pos = (player_pos[0], player_pos[1] - 1)
        elif move == 'R':
            player_pos = (player_pos[0], player_pos[1] + 1)
        elif move == 'U':
            player_pos = (player_pos[0] - 1, player_pos[1])
        elif move == 'D':
            player_pos = (player_pos[0] + 1, player_pos[1])
        
        # Check if player is pushing a box
        if player_pos in boxes:
            box_index = boxes.index(player_pos)
            if move == 'L':
                boxes[box_index] = (boxes[box_index][0], boxes[box_index][1] - 1)
            elif move == 'R':
                boxes[box_index] = (boxes[box_index][0], boxes[box_index][1] + 1)
            elif move == 'U':
                boxes[box_index] = (boxes[box_index][0] - 1, boxes[box_index][1])
            elif move == 'D':
                boxes[box_index] = (boxes[box_index][0] + 1, boxes[box_index][1])
    
    # Check if all boxes are on goals
    if all(box in goals for box in boxes):
        print(moves)
    else:
        print("No solution found")

solve_sokoban()