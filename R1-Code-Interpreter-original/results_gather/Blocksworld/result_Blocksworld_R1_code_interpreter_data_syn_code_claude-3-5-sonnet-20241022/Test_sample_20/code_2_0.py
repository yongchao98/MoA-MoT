def print_solution():
    # Direct sequence of moves to reach the goal state
    moves = [
        "Move H from 2 to 1",  # Move H out of the way
        "Move G from 2 to 3",  # Move G to its final position
        "Move B from 2 to 3",  # Move B temporarily
        "Move E from 2 to 3",  # Move E to its final position
        "Move F from 1 to 2",  # Move F temporarily
        "Move C from 1 to 3",  # Move C to its final position
        "Move A from 3 to 1",  # Move A to its final position
        "Move B from 3 to 1",  # Move B to its final position
        "Move F from 2 to 1"   # Move F to its final position
    ]
    
    print("\n".join(moves))

print_solution()