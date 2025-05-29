def solve_sokoban():
    # Initial positions
    player_pos = (7, 2)
    boxes = [(4, 2), (6, 1), (7, 3)]
    goals = [(2, 2), (3, 1), (7, 4)]
    
    # Moves: R (right), U (up), L (left), D (down)
    moves = "RURRULULDDRR"
    
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
    if all(box in goals for box in boxes):
        print("Solution found:", moves)
    else:
        print("Solution not found")

solve_sokoban()