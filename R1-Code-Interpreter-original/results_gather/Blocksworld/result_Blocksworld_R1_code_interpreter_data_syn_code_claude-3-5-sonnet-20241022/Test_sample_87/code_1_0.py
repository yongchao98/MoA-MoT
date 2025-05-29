def print_solution():
    # Initial state
    stacks = {
        1: ['D', 'B'],
        2: [],
        3: ['C'],
        4: ['A', 'E', 'F']
    }
    
    # Goal state
    goal = {
        1: ['B', 'D'],
        2: ['A', 'C', 'E', 'F'],
        3: []
    }
    
    moves = []
    
    # Step 1: Move B from stack1 to stack2
    moves.append("Move B from 1 to 2")
    stacks[1].pop()  # Remove B
    stacks[2].append('B')
    
    # Step 2: Move D from stack1 to stack3
    moves.append("Move D from 1 to 3")
    stacks[1].pop()  # Remove D
    stacks[3].append('D')
    
    # Step 3: Move F from stack4 to stack1
    moves.append("Move F from 4 to 1")
    stacks[4].pop()  # Remove F
    stacks[1].append('F')
    
    # Step 4: Move E from stack4 to stack3
    moves.append("Move E from 4 to 3")
    stacks[4].pop()  # Remove E
    stacks[3].append('E')
    
    # Step 5: Move A from stack4 to stack2
    moves.append("Move A from 4 to 2")
    stacks[4].pop()  # Remove A
    stacks[2].append('A')
    
    # Step 6: Move E from stack3 to stack2
    moves.append("Move E from 3 to 2")
    stacks[3].pop()  # Remove E
    stacks[2].append('E')
    
    # Step 7: Move F from stack1 to stack2
    moves.append("Move F from 1 to 2")
    stacks[1].pop()  # Remove F
    stacks[2].append('F')
    
    # Step 8: Move D from stack3 to stack1
    moves.append("Move D from 3 to 1")
    stacks[3].pop()  # Remove D
    stacks[1].append('D')
    
    # Step 9: Move C from stack3 to stack2
    moves.append("Move C from 3 to 2")
    stacks[3].pop()  # Remove C
    stacks[2].append('C')
    
    print("<<<" + "\n".join(moves) + ">>>")

print_solution()