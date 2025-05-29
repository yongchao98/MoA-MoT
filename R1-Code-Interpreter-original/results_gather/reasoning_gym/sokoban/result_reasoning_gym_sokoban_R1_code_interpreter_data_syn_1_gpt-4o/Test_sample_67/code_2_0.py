def solve_sokoban():
    # Initial positions
    player_pos = (4, 5)
    boxes = [(2, 6), (3, 3), (4, 4), (5, 4), (6, 5)]
    goals = [(2, 3), (2, 5), (4, 7), (6, 2)]
    
    # Moves: R, D, L, U, R, U, R
    moves = "RDLURUR"
    
    # Simulate the moves
    for move in moves:
        if move == 'R':
            player_pos = (player_pos[0], player_pos[1] + 1)
        elif move == 'L':
            player_pos = (player_pos[0], player_pos[1] - 1)
        elif move == 'U':
            player_pos = (player_pos[0] - 1, player_pos[1])
        elif move == 'D':
            player_pos = (player_pos[0] + 1, player_pos[1])
        
        # Check if player is pushing a box
        if player_pos in boxes:
            box_index = boxes.index(player_pos)
            if move == 'R':
                boxes[box_index] = (boxes[box_index][0], boxes[box_index][1] + 1)
            elif move == 'L':
                boxes[box_index] = (boxes[box_index][0], boxes[box_index][1] - 1)
            elif move == 'U':
                boxes[box_index] = (boxes[box_index][0] - 1, boxes[box_index][1])
            elif move == 'D':
                boxes[box_index] = (boxes[box_index][0] + 1, boxes[box_index][1])
    
    # Check if all boxes are on goals
    solved = all(box in goals for box in boxes)
    return solved

print(solve_sokoban())