def print_solution():
    # Initial state
    stack1 = ['J', 'B']
    stack2 = ['F', 'I', 'E']
    stack3 = ['D', 'A', 'C', 'H', 'G']
    
    # Direct sequence of moves to reach goal state
    moves = []
    
    # First clear the top blocks to make space
    moves.append("Move G from 3 to 1")  # Get G out
    moves.append("Move H from 3 to 2")  # Get H out
    moves.append("Move C from 3 to 1")  # Get C out
    moves.append("Move E from 2 to 3")  # Move E to final position
    moves.append("Move I from 2 to 1")  # Get I out
    moves.append("Move F from 2 to 2")  # Position F
    moves.append("Move D from 3 to 2")  # Position D
    moves.append("Move A from 3 to 3")  # Position A
    moves.append("Move B from 1 to 1")  # Position B
    moves.append("Move C from 1 to 1")  # Position C
    moves.append("Move G from 1 to 1")  # Position G
    moves.append("Move H from 2 to 1")  # Position H
    moves.append("Move I from 1 to 1")  # Position I
    moves.append("Move J from 1 to 1")  # Position J

    print("<<<" + "\n".join(moves) + ">>>")

print_solution()