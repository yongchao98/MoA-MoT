# Simulate the moves to verify the final configuration
def simulate_moves(initial_state, moves):
    state = {k: v[:] for k, v in initial_state.items()}  # Deep copy of the initial state
    for move in moves:
        _, block, _, source, _, destination = move.split()
        block = block.strip()
        source = source.strip()
        destination = destination.strip()
        
        # Perform the move
        state[destination].append(state[source].pop())
    
    return state

# Initial state
initial_state = {
    'stack1': ['E', 'C', 'J'],
    'stack2': ['H', 'D'],
    'stack3': ['G', 'I', 'A', 'B', 'F']
}

# Corrected moves to be performed
moves = [
    "Move F from stack3 to stack1",
    "Move B from stack3 to stack2",
    "Move A from stack3 to stack1",
    "Move G from stack3 to stack1",
    "Move D from stack2 to stack1",
    "Move J from stack1 to stack2"
]

# Simulate the moves
final_state = simulate_moves(initial_state, moves)

# Print the final state
print("Final State:")
for stack, blocks in final_state.items():
    print(f"{stack}: {blocks}")