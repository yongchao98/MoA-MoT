def find_solution():
    # Define initial and goal states
    initial = {
        1: [],
        2: ['B', 'A', 'C'],
        3: ['E', 'D']
    }
    goal = {
        1: ['B', 'E'],
        2: ['C'],
        3: ['A', 'D']
    }
    
    def get_top_block(state, stack):
        return state[stack][-1] if state[stack] else None
        
    def is_valid_move(state, block, from_stack, to_stack):
        # Check if block exists and is on top of from_stack
        if not state[from_stack] or state[from_stack][-1] != block:
            return False
        return True
        
    def apply_move(state, block, from_stack, to_stack):
        new_state = {k: v[:] for k, v in state.items()}
        new_state[from_stack].pop()
        new_state[to_stack].append(block)
        return new_state
        
    def state_to_string(state):
        return str(state)
        
    def is_goal_state(state):
        return all(state[k] == goal[k] for k in state)
    
    def dfs(current_state, moves, visited):
        if is_goal_state(current_state):
            return moves
            
        state_str = state_to_string(current_state)
        if state_str in visited:
            return None
        visited.add(state_str)
        
        # Try all possible moves
        for from_stack in [1, 2, 3]:
            if not current_state[from_stack]:
                continue
            block = current_state[from_stack][-1]
            
            for to_stack in [1, 2, 3]:
                if from_stack != to_stack and is_valid_move(current_state, block, from_stack, to_stack):
                    new_state = apply_move(current_state, block, from_stack, to_stack)
                    new_moves = moves + [f"Move {block} from {from_stack} to {to_stack}"]
                    
                    result = dfs(new_state, new_moves, visited)
                    if result is not None:
                        return result
        
        return None
    
    solution = dfs(initial, [], set())
    if solution:
        print("\n".join(solution))

find_solution()