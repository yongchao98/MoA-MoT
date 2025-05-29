def get_top(stack):
    return stack[-1] if stack else None

def is_valid_move(stacks, from_idx, to_idx):
    if from_idx < 0 or from_idx >= len(stacks) or not stacks[from_idx]:
        return False
    return True

def make_move(stacks, from_idx, to_idx):
    block = stacks[from_idx][-1]
    new_stacks = [stack[:] for stack in stacks]
    new_stacks[from_idx].pop()
    new_stacks[to_idx].append(block)
    return new_stacks

def heuristic(current, goal):
    # Count blocks in wrong positions
    score = 0
    for c_stack, g_stack in zip(current, goal):
        for c_block, g_block in zip(c_stack, g_stack):
            if c_block != g_block:
                score += 1
    return score

def solve_blocks(current, goal, max_depth=15):
    if current == goal:
        return []
    
    stack_count = len(current)
    best_moves = None
    best_score = float('inf')
    
    def dfs(stacks, moves, depth):
        nonlocal best_moves, best_score
        
        if depth == 0 or len(moves) >= max_depth:
            return
        
        if stacks == goal:
            if len(moves) < best_score:
                best_moves = moves[:]
                best_score = len(moves)
            return
        
        current_h = heuristic(stacks, goal)
        if current_h + len(moves) >= best_score:
            return
        
        for from_idx in range(stack_count):
            if not stacks[from_idx]:
                continue
            for to_idx in range(stack_count):
                if from_idx != to_idx:
                    if is_valid_move(stacks, from_idx, to_idx):
                        block = stacks[from_idx][-1]
                        new_stacks = make_move(stacks, from_idx, to_idx)
                        move = f"Move {block} from {from_idx + 1} to {to_idx + 1}"
                        dfs(new_stacks, moves + [move], depth - 1)
    
    dfs(current, [], max_depth)
    return best_moves

# Initial and goal states
initial = [['I', 'A', 'G'], ['D', 'F', 'E'], ['H', 'B', 'C']]
goal = [['B', 'C', 'H'], ['F', 'G', 'I'], ['A', 'D', 'E']]

# Find solution
solution = solve_blocks(initial, goal)

if solution:
    print("<<<")
    print('\n'.join(solution))
    print(">>>")
else:
    print("No solution found")