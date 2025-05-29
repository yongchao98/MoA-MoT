def get_state_string(stacks):
    return '|'.join(','.join(stack) for stack in stacks)

def get_heuristic(stacks):
    goal = [['A', 'C'], ['B', 'D', 'E', 'G'], ['F'], []]
    score = 0
    # Check correct positions
    for i in range(len(stacks)):
        for j in range(len(stacks[i])):
            if i < len(goal) and j < len(goal[i]) and stacks[i][j] == goal[i][j]:
                score += 1
    return -score

def solve_blocksworld():
    initial = [['G'], ['C', 'B'], ['E'], ['A', 'F', 'D']]
    goal = [['A', 'C'], ['B', 'D', 'E', 'G'], ['F'], []]
    
    visited = set()
    states = [(get_heuristic(initial), initial, [])]
    max_iterations = 1000
    iterations = 0
    
    while states and iterations < max_iterations:
        iterations += 1
        _, current, moves = states.pop(0)
        
        if current == goal:
            return moves
            
        state_str = get_state_string(current)
        if state_str in visited:
            continue
        visited.add(state_str)
        
        # Try all possible moves
        for i in range(len(current)):
            if not current[i]:  # Skip empty stacks
                continue
            block = current[i][-1]  # Get top block
            
            for j in range(len(current)):
                if i != j:  # Can't move to same stack
                    # Make move
                    new_state = [stack[:] for stack in current]
                    new_state[i] = new_state[i][:-1]
                    new_state[j] = new_state[j] + [block]
                    
                    new_moves = moves + [f"Move {block} from {i+1} to {j+1}"]
                    h_score = get_heuristic(new_state)
                    states.append((h_score, new_state, new_moves))
        
        # Sort by heuristic score
        states.sort()
        # Keep only the best moves
        states = states[:10]

# Execute solution
solution = solve_blocksworld()
if solution:
    print('\n'.join(solution))