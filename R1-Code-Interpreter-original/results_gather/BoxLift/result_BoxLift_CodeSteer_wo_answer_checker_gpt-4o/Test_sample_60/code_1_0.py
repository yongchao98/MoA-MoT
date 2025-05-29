boxes = [21, 96, 83, 73, 100, 73, 55, 25]
lifters = [62, 79, 54]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Function to check if a combination of lifters can lift a box
def can_lift(box, lifters, used):
    # Try all combinations of lifters
    from itertools import combinations
    for r in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), r):
            if all(not used[i] for i in combo) and sum(lifters[i] for i in combo) >= box:
                return list(combo)
    return None

# Backtracking function
def backtrack(boxes, lifters, step, steps, used):
    if not boxes:
        return True
    if step >= 6:
        return False

    for i, box in enumerate(boxes):
        lifter_indices = can_lift(box, lifters, used)
        if lifter_indices is not None:
            # Mark lifters as used
            for index in lifter_indices:
                used[index] = True
            steps[step].append((box, lifter_indices))
            
            # Debug: Print current step and state
            print(f"Step {step + 1}: Trying to lift box {box} with lifters {lifter_indices}")
            print(f"Remaining boxes: {boxes[:i] + boxes[i+1:]}")
            
            # Recurse to the next step
            if backtrack(boxes[:i] + boxes[i+1:], lifters, step + 1, steps, used):
                return True
            
            # Backtrack: Unmark lifters and remove the step
            steps[step].pop()
            for index in lifter_indices:
                used[index] = False
            
            # Debug: Print backtracking action
            print(f"Backtracking from step {step + 1}, box {box} with lifters {lifter_indices}")

    return False

# Initialize variables
steps = [[] for _ in range(6)]
used = [False] * len(lifters)

# Start backtracking
if backtrack(boxes, lifters, 0, steps, used):
    # Print the steps
    output = "<<<"
    for i, step in enumerate(steps):
        if step:
            output += f"Step {i + 1}: {step}\n"
    output += ">>>"
    print(output)
else:
    print("No solution found within 6 steps.")