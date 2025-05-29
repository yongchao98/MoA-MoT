def get_moves():
    # Initial state representation: each list represents a stack
    initial = [['A', 'E', 'F'], ['B', 'C'], ['D', 'H', 'G'], ['I']]
    goal = [['B', 'C', 'D', 'E'], ['A', 'F', 'G', 'H'], ['I'], []]
    
    # Helper function to get the moves needed
    def solve():
        moves = []
        current = initial[:]
        
        # Step 1: Clear space for final configuration
        # Move F to stack 4
        moves.append("Move F from 1 to 4")
        current[0].pop()  # Remove F
        
        # Move E to stack 3
        moves.append("Move E from 1 to 3")
        current[0].pop()  # Remove E
        
        # Move C to stack 4
        moves.append("Move C from 2 to 4")
        current[1].pop()  # Remove C
        
        # Move B to stack 1
        moves.append("Move B from 2 to 1")
        
        # Move G and H to stack 2
        moves.append("Move G from 3 to 2")
        moves.append("Move H from 3 to 2")
        
        # Move D to stack 1
        moves.append("Move D from 3 to 1")
        
        # Move I to stack 3
        moves.append("Move I from 4 to 3")
        
        # Move C to stack 1
        moves.append("Move C from 4 to 1")
        
        # Move F to stack 2
        moves.append("Move F from 4 to 2")
        
        # Move E to stack 1
        moves.append("Move E from 3 to 1")
        
        # Finally, move A to stack 2
        moves.append("Move A from 1 to 2")
        
        return moves

    # Get the solution
    solution = solve()
    
    # Print in required format
    print("<<<")
    for move in solution:
        print(move)
    print(">>>")

get_moves()