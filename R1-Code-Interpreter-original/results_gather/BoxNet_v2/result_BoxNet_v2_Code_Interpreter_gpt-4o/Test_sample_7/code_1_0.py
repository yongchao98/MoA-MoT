def generate_plan(initial_state, goal_state, adjacency):
    plan = [initial_state.copy()]
    current_state = initial_state.copy()

    while current_state != goal_state:
        for box, current_location in current_state.items():
            if current_location != goal_state[box]:
                for neighbor in adjacency[current_location]:
                    if neighbor == goal_state[box]:
                        current_state[box] = neighbor
                        break
                    elif neighbor not in current_state.values():
                        current_state[box] = neighbor
                        break
        plan.append(current_state.copy())

    return plan

# Define initial and goal states, and adjacency as shown earlier
plan = generate_plan(initial_state, goal_state, adjacency)
print(plan)