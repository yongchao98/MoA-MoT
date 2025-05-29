def print_moves():
    # Initial state
    stack1 = ['D']
    stack2 = ['G', 'C', 'A', 'F']
    stack3 = ['B', 'I', 'E', 'H']
    
    # Store moves
    moves = []
    
    # Step 1: Move F from stack2 to clear A
    moves.append("Move F from 2 to 3")
    
    # Step 2: Move A to its final position in stack1
    moves.append("Move A from 2 to 1")
    
    # Step 3: Move F to its final position in stack1
    moves.append("Move F from 3 to 1")
    
    # Step 4: Move H to its final position in stack1
    moves.append("Move H from 3 to 1")
    
    # Step 5: Move E to clear I
    moves.append("Move E from 3 to 2")
    
    # Step 6: Move I to its final position in stack2
    moves.append("Move I from 3 to 2")
    
    # Step 7: Move G to its final position in stack2
    moves.append("Move G from 2 to 2")
    
    # Step 8: Move C to its final position in stack2
    moves.append("Move C from 2 to 2")
    
    # Step 9: Move D to its final position in stack3
    moves.append("Move D from 1 to 3")
    
    # Step 10: Move E to its final position in stack3
    moves.append("Move E from 2 to 3")
    
    # Print the moves in required format
    print("<<<")
    print("\n".join(moves))
    print(">>>")

print_moves()