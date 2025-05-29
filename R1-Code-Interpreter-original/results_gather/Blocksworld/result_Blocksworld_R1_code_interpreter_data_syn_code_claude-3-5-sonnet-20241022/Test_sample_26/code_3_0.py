def print_solution():
    moves = []
    # First get C to stack 1
    moves.append("Move C from 1 to 3")  # Save C temporarily
    moves.append("Move J from 1 to 2")  # Move J towards goal position
    moves.append("Move K from 1 to 3")  # Move K out of the way
    moves.append("Move B from 1 to 2")  # Move B towards goal position
    
    # Get A to stack 3 (its goal position)
    moves.append("Move A from 3 to 3")  # Keep A in stack 3
    
    # Build stack 1 (C D F I K)
    moves.append("Move C from 3 to 1")  # Start building stack 1 from bottom
    moves.append("Move D from 3 to 1")  # Add D
    moves.append("Move F from 2 to 1")  # Add F
    moves.append("Move I from 3 to 1")  # Add I
    moves.append("Move K from 3 to 1")  # Add K
    
    # Stack 2 is already close to goal (B H J)
    moves.append("Move H from 2 to 2")  # Keep H in position
    moves.append("Move J from 2 to 2")  # Keep J in position
    
    # Stack 3 final adjustments for (A E G)
    moves.append("Move E from 3 to 3")  # Keep E in position
    moves.append("Move G from 3 to 3")  # Keep G in position
    
    print('\n'.join(moves))

print_solution()