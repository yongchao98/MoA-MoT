def print_solution():
    # Initial state
    stacks = [['A', 'E'], ['D'], ['C', 'B']]
    moves = []
    
    # Step 1: Move B from stack 3 to stack 1
    moves.append("Move B from 3 to 1")
    
    # Step 2: Move C from stack 3 to stack 2
    moves.append("Move C from 3 to 2")
    
    # Step 3: Move D from stack 2 to stack 3
    moves.append("Move D from 2 to 3")
    
    # Print solution
    print("<<<" + "\n".join(moves) + ">>>")

print_solution()