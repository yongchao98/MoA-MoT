def print_solution():
    # Initial state
    initial = {
        1: ['K', 'I', 'G'],
        2: ['A', 'F', 'B'],
        3: ['J', 'E'],
        4: ['C', 'H', 'D']
    }
    
    # Direct sequence of moves to reach goal state
    moves = [
        # Clear stack 4 first
        "Move D from 4 to 2",
        "Move H from 4 to 1",
        "Move C from 4 to 3",
        
        # Arrange stack 1
        "Move G from 1 to 3",
        "Move I from 1 to 2",
        "Move K from 1 to 3",
        
        # Move blocks to their final positions
        "Move B from 2 to 1",
        "Move F from 2 to 1",
        "Move E from 3 to 1",
        "Move J from 3 to 3"
    ]
    
    print("\n".join(moves))

print_solution()