def solve_sokoban():
    # Initial positions
    player_pos = (5, 2)
    boxes = [(1, 2), (2, 4), (4, 2)]
    goals = [(3, 2), (5, 3)]
    
    # Moves: L, R, U, D
    moves = []
    
    # Move player to (4, 2) and push box to (3, 2)
    moves.extend(['U', 'U', 'R', 'R', 'D', 'D', 'L', 'U'])
    
    # Move player to (2, 4) and push box to (5, 3)
    moves.extend(['D', 'R', 'R', 'U', 'U', 'L', 'D', 'D', 'R'])
    
    # Move player to (1, 2) and push box to (5, 3)
    moves.extend(['U', 'U', 'L', 'L', 'D', 'D', 'R', 'R', 'D', 'D'])
    
    # Convert moves to a string
    solution = ''.join(moves)
    print(solution)

solve_sokoban()