def print_solution():
    # Initial state
    stack1 = ['D']
    stack2 = ['A', 'J']
    stack3 = ['H', 'F', 'C', 'E', 'I']
    stack4 = ['B', 'G']
    
    moves = []
    
    # Move I from stack3 to stack1
    moves.append("Move I from 3 to 1")
    stack3.pop()  # Remove I
    stack1.append('I')
    
    # Move E from stack3 to stack4
    moves.append("Move E from 3 to 4")
    stack3.pop()  # Remove E
    stack4.append('E')
    
    # Move C from stack3 to stack4
    moves.append("Move C from 3 to 4")
    stack3.pop()  # Remove C
    stack4.append('C')
    
    # Move F from stack3 to stack1
    moves.append("Move F from 3 to 1")
    stack3.pop()  # Remove F
    stack1.append('F')
    
    # Move J from stack2 to stack3
    moves.append("Move J from 2 to 3")
    stack2.pop()  # Remove J
    stack3.append('J')
    
    # Move G from stack4 to stack3
    moves.append("Move G from 4 to 3")
    stack4.pop()  # Remove G
    stack3.append('G')
    
    # Move C from stack4 to stack2
    moves.append("Move C from 4 to 2")
    stack4.pop()  # Remove C
    stack2.append('C')
    
    # Move D from stack1 to stack2
    moves.append("Move D from 1 to 2")
    stack1.pop()  # Remove D
    stack2.append('D')
    
    # Move G from stack3 to stack2
    moves.append("Move G from 3 to 2")
    stack3.pop()  # Remove G
    stack2.append('G')
    
    # Move E from stack4 to stack3
    moves.append("Move E from 4 to 3")
    stack4.pop()  # Remove E
    stack3.append('E')
    
    # Move H to be in correct position in stack3
    moves.append("Move H from 3 to 3")
    
    print("<<<")
    print("\n".join(moves))
    print(">>>")

print_solution()