from collections import deque
import time

def blocks_world_bfs(initial_state, goal_state, timeout=10):
    start_time = time.time()
    queue = deque([(initial_state, [])])
    visited = set()
    
    while queue:
        # Check for timeout
        if time.time() - start_time > timeout:
            print("Timeout: Solution not found within the time limit.")
            return
        
        current_state, moves = queue.popleft()
        
        if current_state == goal_state:
            print("<<<" + "\n".join(moves) + ">>>")
            return
        
        visited.add(current_state)
        
        for i, stack in enumerate(current_state):
            if stack:
                block = stack[-1]
                for j in range(len(current_state)):
                    if i != j:
                        new_state = list(map(list, current_state))
                        new_state[i].pop()
                        new_state[j].append(block)
                        new_state_tuple = tuple(map(tuple, new_state))
                        
                        if new_state_tuple not in visited:
                            new_moves = moves + [f"Move {block} from stack{i+1} to stack{j+1}"]
                            queue.append((new_state_tuple, new_moves))

# Initial and goal states
initial_state = (('B', 'K', 'J', 'C'), ('F', 'H'), ('D', 'E', 'G', 'I', 'A'))
goal_state = (('C', 'D', 'F', 'I', 'K'), ('B', 'H', 'J'), ('A', 'E', 'G'))

# Find the sequence of moves
blocks_world_bfs(initial_state, goal_state)