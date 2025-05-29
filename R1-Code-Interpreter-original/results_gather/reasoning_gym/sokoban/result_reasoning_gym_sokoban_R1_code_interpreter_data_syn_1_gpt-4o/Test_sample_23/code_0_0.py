def solve_sokoban():
    # Initial positions
    player_pos = (1, 4)
    boxes = [(1, 3), (2, 4), (3, 4)]
    goals = [(1, 2), (3, 5), (4, 5)]
    
    # Moves: L, R, U, D
    moves = []
    
    # Move 1: Push box from (1, 3) to (1, 2)
    moves.extend(['L', 'L'])
    
    # Move 2: Move player to (2, 4) and push box to (3, 5)
    moves.extend(['R', 'R', 'D', 'D', 'L', 'U', 'R', 'R', 'D'])
    
    # Move 3: Move player to (3, 4) and push box to (4, 5)
    moves.extend(['U', 'L', 'D', 'D', 'R'])
    
    # Return the solution as a string
    return ''.join(moves)

# Print the solution
print(solve_sokoban())