def get_moves():
    # Initial state representation
    s1 = ['F', 'A']
    s2 = ['K', 'D', 'E', 'B', 'I', 'C']
    s3 = ['H', 'J', 'G']
    
    # Goal state
    g1 = ['G']
    g2 = ['B', 'D']
    g3 = ['A', 'C', 'E', 'F', 'H', 'I', 'J', 'K']
    
    moves = []
    
    # Direct approach based on goal state analysis
    # First clear path for G to reach stack1
    moves.append(f"Move C from 2 to 3")  # Move C to its final position
    moves.append(f"Move I from 2 to 3")  # Move I to its final position
    moves.append(f"Move B from 2 to 2")  # Keep B for stack2
    moves.append(f"Move E from 2 to 3")  # Move E to its final position
    moves.append(f"Move D from 2 to 2")  # Keep D for stack2
    moves.append(f"Move K from 2 to 3")  # Move K to its final position
    
    # Move blocks from stack1 to their final positions
    moves.append(f"Move A from 1 to 3")  # Move A to its final position
    moves.append(f"Move F from 1 to 3")  # Move F to its final position
    
    # Move G to its final position
    moves.append(f"Move G from 3 to 1")  # G goes to stack1
    
    # Final adjustments for stack3
    moves.append(f"Move J from 3 to 3")  # Ensure J is in position
    moves.append(f"Move H from 3 to 3")  # Ensure H is in position
    
    print("<<<" + "\n".join(moves) + ">>>")

get_moves()