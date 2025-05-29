def solve_sokoban():
    # Initial positions
    player_pos = (1, 1)
    boxes = [(1, 2), (2, 1), (2, 2)]
    goals = [(3, 1), (3, 2)]
    
    # Moves: U, D, L, R
    moves = []
    
    # Move 1: Push box from (2, 1) to (3, 1)
    moves.extend(['D', 'R', 'R', 'U', 'L', 'D'])
    
    # Move 2: Push box from (2, 2) to (3, 2)
    moves.extend(['U', 'R', 'D', 'L', 'D'])
    
    # Convert moves to a string
    solution = ''.join(moves)
    print(solution)

solve_sokoban()