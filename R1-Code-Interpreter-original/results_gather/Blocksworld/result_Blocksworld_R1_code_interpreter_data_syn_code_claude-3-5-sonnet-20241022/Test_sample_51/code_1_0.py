def get_moves():
    # Initial state representation
    stacks = {
        1: ['G', 'A', 'H'],
        2: ['E'],
        3: ['I', 'C', 'B'],
        4: ['F', 'D']
    }
    
    # Goal state for checking
    goal = {
        1: ['A', 'B', 'C', 'D', 'G', 'H', 'I'],
        2: ['E'],
        3: ['F'],
        4: []
    }
    
    moves = []
    
    # Step 1: Clear space for final configuration
    # Move H from top of stack 1
    moves.append(f"Move H from 1 to 4")
    stacks[4].append(stacks[1].pop())
    
    # Move B from stack 3
    moves.append(f"Move B from 3 to 1")
    stacks[1].append(stacks[3].pop())
    
    # Move D from stack 4
    moves.append(f"Move D from 4 to 2")
    stacks[2].append(stacks[4].pop())
    
    # Move F to its final position
    moves.append(f"Move F from 4 to 3")
    stacks[3].append(stacks[4].pop())
    
    # Build up stack 1
    moves.append(f"Move C from 3 to 1")
    stacks[1].append(stacks[3].pop())
    
    moves.append(f"Move I from 3 to 1")
    stacks[1].append(stacks[3].pop())
    
    moves.append(f"Move D from 2 to 1")
    stacks[1].append(stacks[2].pop())
    
    # Move G temporarily
    moves.append(f"Move G from 1 to 4")
    stacks[4].append(stacks[1].pop())
    
    # Move H to final position
    moves.append(f"Move H from 4 to 1")
    stacks[1].append(stacks[4].pop())
    
    # Move G to final position
    moves.append(f"Move G from 4 to 1")
    stacks[1].append(stacks[4].pop())
    
    print("<<<" + "\n".join(moves) + ">>>")

get_moves()