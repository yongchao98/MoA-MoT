from collections import deque

def blocksworld_bfs(initial_state, goal_state):
    # Queue for BFS: stores (current_state, path_to_state)
    queue = deque([(initial_state, [])])
    visited = set()
    
    while queue:
        current_state, path = queue.popleft()
        
        # Check if the current state is the goal state
        if current_state == goal_state:
            return path
        
        # Mark the current state as visited
        visited.add(current_state)
        
        # Generate all possible moves
        for i, stack in enumerate(current_state):
            if stack:  # If the stack is not empty
                block = stack[-1]  # Top block
                for j in range(len(current_state)):
                    if i != j:  # Move to a different stack
                        new_state = list(map(list, current_state))
                        new_state[i].pop()  # Remove block from current stack
                        new_state[j].append(block)  # Add block to new stack
                        new_state_tuple = tuple(map(tuple, new_state))
                        
                        if new_state_tuple not in visited:
                            move = f"Move {block} from {i+1} to {j+1}"
                            queue.append((new_state_tuple, path + [move]))
    
    return None  # No solution found

# Initial and goal states
initial_state = ((), ('B', 'E'), ('C', 'F', 'D'), ('A',))
goal_state = (('A', 'B', 'F'), ('C',), ('D', 'E'), ())

# Find the sequence of moves
moves = blocksworld_bfs(initial_state, goal_state)

# Output the moves
if moves:
    print("<<<" + "\n".join(moves) + ">>>")
else:
    print("No solution found.")