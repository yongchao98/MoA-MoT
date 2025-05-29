def print_solution():
    # Initial state
    stack1 = ['D']
    stack2 = ['A', 'J']
    stack3 = ['H', 'F', 'C', 'E', 'I']
    stack4 = ['B', 'G']
    
    moves = []
    
    # Move I from stack3 to stack1
    moves.append("Move I from 3 to 1")
    
    # Move E from stack3 to stack4
    moves.append("Move E from 3 to 4")
    
    # Move C from stack3 to stack4
    moves.append("Move C from 3 to 4")
    
    # Move F from stack3 to stack1
    moves.append("Move F from 3 to 1")
    
    # Move J from stack2 to stack3
    moves.append("Move J from 2 to 3")
    
    # Move G from stack4 to stack2
    moves.append("Move G from 4 to 2")
    
    # Move B from stack4 to stack2
    moves.append("Move B from 4 to 2")
    
    # Move C from stack4 to stack2
    moves.append("Move C from 4 to 2")
    
    # Move D from stack1 to stack2
    moves.append("Move D from 1 to 2")
    
    # Move E from stack4 to stack3
    moves.append("Move E from 3 to 3")
    
    # Move H to be in correct position in stack3
    moves.append("Move H from 3 to 3")
    
    print("<<<")
    print("\n".join(moves))
    print(">>>")

print_solution()