def get_top_block(stack):
    return stack[-1] if stack else None

def is_valid_move(stacks, from_stack, to_stack):
    if not stacks[from_stack]:  # source stack is empty
        return False
    return True

def make_move(state, move_from, move_to):
    new_state = [list(stack) for stack in state]
    block = new_state[move_from].pop()
    new_state[move_to].append(block)
    return new_state

def solve_blocks():
    # Initial and goal states (0-based indexing for stacks)
    initial = [['B', 'G'], ['F', 'E', 'A', 'C'], ['H', 'D']]
    goal = [['A', 'B', 'F'], ['D', 'E'], ['C', 'G', 'H']]
    
    def dfs(current_state, moves, depth):
        if depth > 15:  # Limit depth to prevent excessive searching
            return None
        
        if current_state == goal:
            return moves
            
        for from_stack in range(3):
            for to_stack in range(3):
                if from_stack != to_stack and current_state[from_stack]:
                    # Check if move is valid
                    if is_valid_move(current_state, from_stack, to_stack):
                        # Make the move
                        block = current_state[from_stack][-1]
                        new_state = make_move(current_state, from_stack, to_stack)
                        new_move = f"Move {block} from {from_stack + 1} to {to_stack + 1}"
                        
                        # Recursive call
                        result = dfs(new_state, moves + [new_move], depth + 1)
                        if result is not None:
                            return result
        return None

    solution = dfs(initial, [], 0)
    if solution:
        print('\n'.join(solution))

solve_blocks()