def get_top_block(stack):
    return stack[-1] if stack else None

def print_moves(initial_stacks, goal_stacks, max_depth=20):
    def is_valid_move(stacks, from_idx, to_idx):
        return bool(stacks[from_idx])  # Source stack must not be empty

    def make_move(stacks, from_idx, to_idx):
        new_stacks = [stack[:] for stack in stacks]
        block = new_stacks[from_idx].pop()
        new_stacks[to_idx].append(block)
        return new_stacks

    def dfs(current, moves, depth):
        if depth > max_depth:
            return None
        
        if current == goal_stacks:
            return moves
        
        for i in range(3):  # From stack
            for j in range(3):  # To stack
                if i != j and is_valid_move(current, i, j):
                    block = get_top_block(current[i])
                    new_stacks = make_move(current, i, j)
                    new_move = f"Move {block} from {i+1} to {j+1}"
                    
                    result = dfs(new_stacks, moves + [new_move], depth + 1)
                    if result:
                        return result
        return None

    # Initial and goal states
    initial = [['H', 'D'], 
              ['E', 'J', 'A', 'C', 'F'],
              ['B', 'I', 'G', 'K']]
    
    goal = [['B', 'E', 'I', 'J', 'K'],
            ['A', 'C', 'D', 'F'],
            ['G', 'H']]

    solution = dfs(initial, [], 0)
    
    if solution:
        print("<<<")
        print('\n'.join(solution))
        print(">>>")
    else:
        print("No solution found within depth limit")

# Run the solution
print_moves([['H', 'D'], 
            ['E', 'J', 'A', 'C', 'F'],
            ['B', 'I', 'G', 'K']],
           [['B', 'E', 'I', 'J', 'K'],
            ['A', 'C', 'D', 'F'],
            ['G', 'H']])