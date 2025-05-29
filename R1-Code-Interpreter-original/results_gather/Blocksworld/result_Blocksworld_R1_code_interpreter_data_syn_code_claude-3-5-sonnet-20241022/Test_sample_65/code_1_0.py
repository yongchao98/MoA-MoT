def get_moves():
    # Initial state
    stacks = [
        ['B', 'G', 'I'],    # stack 1
        ['F', 'H', 'D'],    # stack 2
        ['E', 'C'],         # stack 3
        ['A']               # stack 4
    ]
    
    # Goal state
    goal = [
        ['B', 'D', 'E', 'I'],  # stack 1
        ['F', 'G', 'H'],       # stack 2
        ['A', 'C'],            # stack 3
        []                     # stack 4
    ]
    
    moves = []
    
    # Direct solution based on analysis
    # Move I to make space
    moves.append(f"Move I from 1 to 4")
    stacks[0].pop()  # Remove I
    stacks[3].append('I')
    
    # Move C to prepare for A
    moves.append(f"Move C from 3 to 1")
    stacks[2].pop()  # Remove C
    stacks[0].append('C')
    
    # Move A to its final position
    moves.append(f"Move A from 4 to 3")
    stacks[3].pop()  # Remove A
    stacks[2].append('A')
    
    # Move D to its position
    moves.append(f"Move D from 2 to 1")
    stacks[1].pop()  # Remove D
    stacks[0].append('D')
    
    # Move H temporarily
    moves.append(f"Move H from 2 to 4")
    stacks[1].pop()  # Remove H
    stacks[3].append('H')
    
    # Move G to its final position
    moves.append(f"Move G from 1 to 2")
    stacks[0].pop()  # Remove G
    stacks[1].append('G')
    
    # Move H to its final position
    moves.append(f"Move H from 4 to 2")
    stacks[3].pop()  # Remove H
    stacks[1].append('H')
    
    # Move E to its position
    moves.append(f"Move E from 3 to 1")
    stacks[2].pop()  # Remove E
    stacks[0].append('E')
    
    # Move I back to its final position
    moves.append(f"Move I from 4 to 1")
    stacks[3].pop()  # Remove I
    stacks[0].append('I')
    
    print('\n'.join(moves))

get_moves()