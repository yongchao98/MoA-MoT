def get_solution():
    # Initial state analysis:
    # stack1: J
    # stack2: D H F I B A G (G is on top)
    # stack3: E C (C is on top)
    
    # Goal state:
    # stack1: G I
    # stack2: A B C D E F H
    # stack3: J

    moves = []
    
    # Step 1: Move G to stack1 (it's already on top of stack2)
    moves.append("Move G from 2 to 1")
    
    # Step 2: Move A from stack2 to stack3 (clearing the way to I)
    moves.append("Move A from 2 to 3")
    
    # Step 3: Move B from stack2 to stack3
    moves.append("Move B from 2 to 3")
    
    # Step 4: Move I from stack2 to stack1
    moves.append("Move I from 2 to 1")
    
    # Step 5: Move remaining blocks from stack2 to stack3 temporarily
    moves.append("Move F from 2 to 3")
    moves.append("Move H from 2 to 3")
    moves.append("Move D from 2 to 3")
    
    # Step 6: Move J from stack1 to stack3
    moves.append("Move J from 1 to 3")
    
    # Step 7: Build stack2 in correct order
    moves.append("Move D from 3 to 2")
    moves.append("Move H from 3 to 2")
    moves.append("Move F from 3 to 2")
    moves.append("Move E from 3 to 2")
    moves.append("Move C from 3 to 2")
    moves.append("Move B from 3 to 2")
    moves.append("Move A from 3 to 2")
    
    print("<<<" + "\n".join(moves) + ">>>")

get_solution()