boxes = [99, 159, 55, 86, 215, 56, 271, 167, 212, 87, 78, 211, 246, 153, 227, 287]
lifters = [143, 50, 70, 79, 136]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Function to find a combination of lifters to lift a box
def find_combination(box, lifters, used):
    if box == 0:
        return []
    for i, lifter in enumerate(lifters):
        if not used[i] and lifter <= box:
            used[i] = True
            result = find_combination(box - lifter, lifters, used)
            if result is not None:
                return [i] + result
            used[i] = False
    return None

# Backtracking function to assign lifters to boxes
def assign_lifters(boxes, lifters, step, max_steps, steps):
    if not boxes:
        return True
    if step >= max_steps:
        return False
    
    used = [False] * len(lifters)
    current_step = []
    
    for box in boxes[:]:
        combination = find_combination(box, lifters, used)
        if combination is not None:
            current_step.append((box, combination))
            boxes.remove(box)
    
    if current_step:
        steps.append(current_step)
        if assign_lifters(boxes, lifters, step + 1, max_steps, steps):
            return True
        steps.pop()
        boxes.extend([box for box, _ in current_step])
    
    return False

# Main function to solve the problem
def solve_box_lift(boxes, lifters, max_steps):
    steps = []
    if assign_lifters(boxes, lifters, 0, max_steps, steps):
        output = []
        for i, step in enumerate(steps):
            output.append(f"Step {i + 1}: {step}")
        print("<<<" + "\n".join(output) + ">>>")
    else:
        print("No solution found within the given steps.")

solve_box_lift(boxes, lifters, 7)