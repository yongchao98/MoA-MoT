def get_state_distance(state, goal):
    # Simple heuristic: count blocks in wrong positions
    distance = 0
    for i in range(len(state)):
        for j in range(len(state[i])):
            if j < len(goal[i]) and state[i][j] != goal[i][j]:
                distance += 1
            elif j >= len(goal[i]):
                distance += 1
    return distance

def solve_blocksworld():
    initial_state = [['D', 'I'], ['H', 'A', 'C'], ['F', 'B', 'E'], ['G']]
    goal_state = [['A', 'B', 'C', 'D', 'F'], ['E', 'H', 'I'], ['G'], []]
    
    # Priority queue based on heuristic distance
    states = [(get_state_distance(initial_state, goal_state), initial_state, [])]
    visited = {str(initial_state)}
    max_iterations = 1000  # Limit search space
    
    for _ in range(max_iterations):
        if not states:
            break
            
        _, current_state, moves = states.pop(0)
        
        if current_state == goal_state:
            return moves
            
        # Try moves that might be helpful
        for from_stack in range(4):
            if not current_state[from_stack]:
                continue
                
            block = current_state[from_stack][-1]
            # Find where this block should be in goal state
            target_stack = None
            for i, stack in enumerate(goal_state):
                if block in stack:
                    target_stack = i
                    break
            
            # Prioritize moves towards goal position
            possible_moves = [target_stack] if target_stack is not None else range(4)
            for to_stack in possible_moves:
                if from_stack != to_stack:
                    new_state = [stack[:] for stack in current_state]
                    new_state[from_stack] = new_state[from_stack][:-1]
                    new_state[to_stack] = new_state[to_stack] + [block]
                    
                    state_str = str(new_state)
                    if state_str not in visited:
                        visited.add(state_str)
                        new_moves = moves + [f"Move {block} from {from_stack + 1} to {to_stack + 1}"]
                        distance = get_state_distance(new_state, goal_state)
                        # Insert maintaining sort by distance
                        insert_pos = 0
                        while insert_pos < len(states) and states[insert_pos][0] < distance:
                            insert_pos += 1
                        states.insert(insert_pos, (distance, new_state, new_moves))
    
    return None

solution = solve_blocksworld()
if solution:
    print('\n'.join(solution))