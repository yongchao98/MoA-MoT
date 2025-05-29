def solve_sokoban():
    # Initial positions
    player_pos = (1, 4)
    boxes = [(2, 2), (2, 3), (3, 5)]
    goals = [(1, 2), (4, 4)]
    
    # Moves: L, R, U, D
    moves = []
    
    # Move player to push the first box at (2, 2) to the goal at (1, 2)
    moves.extend(['U', 'L', 'L', 'D', 'R', 'U'])
    
    # Now, the player is at (1, 3), and the box is on the goal at (1, 2)
    # Move player to push the second box at (2, 3) to the goal at (4, 4)
    moves.extend(['D', 'R', 'R', 'D', 'D', 'L', 'U', 'L', 'D', 'R', 'R', 'U', 'L', 'D'])
    
    # Now, the player is at (4, 5), and the second box is on the goal at (4, 4)
    # Move player to push the third box at (3, 5) to the goal at (4, 4)
    moves.extend(['U', 'L', 'L', 'D', 'R', 'R', 'U', 'L', 'D'])
    
    # Convert moves to a string
    return ''.join(moves)

# Print the solution
print(solve_sokoban())