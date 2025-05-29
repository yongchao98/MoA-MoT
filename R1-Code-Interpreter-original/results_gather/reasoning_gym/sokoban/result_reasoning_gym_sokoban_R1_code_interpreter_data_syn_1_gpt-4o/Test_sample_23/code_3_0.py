def sokoban_solver():
    # Initial positions
    player_pos = (1, 4)
    boxes = [(1, 3), (2, 4), (3, 5)]
    goals = [(1, 2), (3, 4), (4, 5)]
    
    # Moves: L, D, D
    moves = "LDD"
    
    # Simulate the moves
    for move in moves:
        if move == 'R':
            # Move player right
            player_pos = (player_pos[0], player_pos[1] + 1)
            # Check if player is pushing a box
            if player_pos in boxes:
                # Move the box
                box_index = boxes.index(player_pos)
                boxes[box_index] = (boxes[box_index][0], boxes[box_index][1] + 1)
        elif move == 'D':
            # Move player down
            player_pos = (player_pos[0] + 1, player_pos[1])
            # Check if player is pushing a box
            if player_pos in boxes:
                # Move the box
                box_index = boxes.index(player_pos)
                boxes[box_index] = (boxes[box_index][0] + 1, boxes[box_index][1])
        elif move == 'L':
            # Move player left
            player_pos = (player_pos[0], player_pos[1] - 1)
            # Check if player is pushing a box
            if player_pos in boxes:
                # Move the box
                box_index = boxes.index(player_pos)
                boxes[box_index] = (boxes[box_index][0], boxes[box_index][1] - 1)
        elif move == 'U':
            # Move player up
            player_pos = (player_pos[0] - 1, player_pos[1])
            # Check if player is pushing a box
            if player_pos in boxes:
                # Move the box
                box_index = boxes.index(player_pos)
                boxes[box_index] = (boxes[box_index][0] - 1, boxes[box_index][1])
    
    # Check if all boxes are on goals
    solved = all(box in goals for box in boxes)
    return solved

print(sokoban_solver())