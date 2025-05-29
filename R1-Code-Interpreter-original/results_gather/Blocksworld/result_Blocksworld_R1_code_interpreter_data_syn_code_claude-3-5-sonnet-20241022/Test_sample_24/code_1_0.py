def move_block(stacks, moves, depth=0, max_depth=15):
    if depth >= max_depth:
        return False
    
    # Check if we reached the goal state
    if (stacks[0] == ['E', 'H'] and 
        stacks[1] == ['A', 'G'] and 
        stacks[2] == ['B', 'C', 'D', 'F', 'I', 'J']):
        return True
    
    # Try all possible moves
    for from_stack in range(3):
        if not stacks[from_stack]:  # Skip empty stacks
            continue
        
        block = stacks[from_stack][-1]  # Get top block
        
        for to_stack in range(3):
            if from_stack == to_stack:
                continue
                
            # Make move
            stacks[from_stack].pop()
            stacks[to_stack].append(block)
            moves.append(f"Move {block} from {from_stack + 1} to {to_stack + 1}")
            
            # Recursively try next moves
            if move_block(stacks, moves, depth + 1, max_depth):
                return True
            
            # Undo move if it didn't lead to solution
            stacks[to_stack].pop()
            stacks[from_stack].append(block)
            moves.pop()
    
    return False

# Initial state
stacks = [
    ['I', 'D', 'B'],
    ['G', 'J', 'A', 'C', 'E'],
    ['H', 'F']
]

moves = []
if move_block(stacks, moves):
    print("<<<")
    print("\n".join(moves))
    print(">>>")
else:
    print("No solution found")