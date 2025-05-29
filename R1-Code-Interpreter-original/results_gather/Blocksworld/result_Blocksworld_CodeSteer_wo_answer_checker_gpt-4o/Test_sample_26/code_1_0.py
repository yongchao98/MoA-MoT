import heapq

def heuristic(state, goal_state):
    # Calculate the heuristic cost as the number of misplaced blocks
    cost = 0
    for stack_index, stack in enumerate(state):
        for block_index, block in enumerate(stack):
            if block not in goal_state[stack_index] or goal_state[stack_index].index(block) != block_index:
                cost += 1
    return cost

def blocks_world_astar(initial_state, goal_state):
    # Priority queue for A* search: stores tuples of (estimated_cost, current_state, moves)
    priority_queue = []
    heapq.heappush(priority_queue, (0, initial_state, []))
    visited = set()

    while priority_queue:
        estimated_cost, current_state, moves = heapq.heappop(priority_queue)

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
                            cost = len(new_moves) + heuristic(new_state_tuple, goal_state)
                            heapq.heappush(priority_queue, (cost, new_state_tuple, new_moves))

# Initial and goal states
initial_state = (('B', 'K', 'J', 'C'), ('F', 'H'), ('D', 'E', 'G', 'I', 'A'))
goal_state = (('C', 'D', 'F', 'I', 'K'), ('B', 'H', 'J'), ('A', 'E', 'G'))

# Find the sequence of moves
blocks_world_astar(initial_state, goal_state)