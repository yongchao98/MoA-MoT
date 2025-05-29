def is_goal_state(stacks, goal_stacks):
    return stacks == goal_stacks

def get_possible_moves(stacks):
    moves = []
    for i, stack in enumerate(stacks):
        if stack:  # If the stack is not empty
            block = stack[-1]
            for j, target_stack in enumerate(stacks):
                if i != j:  # Don't move to the same stack
                    moves.append((block, i, j))
    return moves

def apply_move(stacks, move):
    block, source, destination = move
    stacks[source].pop()
    stacks[destination].append(block)

def backtrack(stacks, goal_stacks, moves_sequence, visited_states):
    if is_goal_state(stacks, goal_stacks):
        return True

    state_tuple = tuple(tuple(stack) for stack in stacks)
    if state_tuple in visited_states:
        return False

    visited_states.add(state_tuple)

    for move in get_possible_moves(stacks):
        apply_move(stacks, move)
        moves_sequence.append(move)
        
        if backtrack(stacks, goal_stacks, moves_sequence, visited_states):
            return True
        
        # Undo the move
        moves_sequence.pop()
        apply_move(stacks, (move[0], move[2], move[1]))

    visited_states.remove(state_tuple)
    return False

def solve_blocksworld(initial_stacks, goal_stacks):
    moves_sequence = []
    visited_states = set()
    if backtrack(initial_stacks, goal_stacks, moves_sequence, visited_states):
        return moves_sequence
    else:
        return None

# Initial and goal states
initial_stacks = [['H', 'J', 'A', 'B'], ['F', 'I', 'K', 'C'], [], ['G', 'E', 'D']]
goal_stacks = [['B', 'C', 'D', 'E', 'I', 'K'], ['A', 'H', 'J'], ['F', 'G'], []]

# Solve the problem
solution = solve_blocksworld(initial_stacks, goal_stacks)

# Print the solution
if solution:
    print("<<<")
    for move in solution:
        print(f"Move {move[0]} from {move[1] + 1} to {move[2] + 1}")
    print(">>>")
else:
    print("No solution found.")