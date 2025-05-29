from collections import deque

def is_goal_state(stacks, goal):
    return stacks == goal

def get_possible_moves(stacks):
    moves = []
    for i, stack in enumerate(stacks):
        if stack:  # If the stack is not empty
            block = stack[-1]
            for j in range(len(stacks)):
                if i != j:  # Move to a different stack
                    new_stacks = [list(s) for s in stacks]  # Deep copy
                    new_stacks[i].pop()
                    new_stacks[j].append(block)
                    moves.append((block, i + 1, j + 1, new_stacks))
    return moves

def bfs(initial_state, goal_state):
    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(tuple(tuple(s) for s in initial_state))
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_goal_state(current_state, goal_state):
            return path
        
        for block, src, dest, new_stacks in get_possible_moves(current_state):
            state_tuple = tuple(tuple(s) for s in new_stacks)
            if state_tuple not in visited:
                visited.add(state_tuple)
                queue.append((new_stacks, path + [f"Move {block} from stack{src} to stack{dest}"]))
    
    return None

initial_state = [['H', 'G', 'A', 'C'], ['E'], ['F'], ['B', 'D']]
goal_state = [['A', 'B', 'E', 'G'], ['C', 'F'], ['D', 'H']]

solution = bfs(initial_state, goal_state)
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")