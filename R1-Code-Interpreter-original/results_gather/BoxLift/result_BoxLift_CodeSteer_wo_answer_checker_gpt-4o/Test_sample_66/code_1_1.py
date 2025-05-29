def can_lift(box, lifters, used):
    lifter_indices = []
    remaining_weight = box
    for i, lifter in enumerate(lifters):
        if not used[i] and lifter <= remaining_weight:
            lifter_indices.append(i)
            remaining_weight -= lifter
            if remaining_weight <= 0:
                return lifter_indices
    return None

def backtrack(boxes, lifters, step, max_steps, used, current_solution, best_solution):
    if step >= max_steps:
        return False
    if not boxes:
        if len(current_solution) < len(best_solution[0]):
            best_solution[0] = current_solution[:]
        return True
    
    for i, box in enumerate(boxes):
        lifter_indices = can_lift(box, lifters, used)
        if lifter_indices is not None:
            for idx in lifter_indices:
                used[idx] = True
            current_solution.append((box, lifter_indices))
            remaining_boxes = boxes[:i] + boxes[i+1:]
            if backtrack(remaining_boxes, lifters, step + 1, max_steps, used, current_solution, best_solution):
                return True
            current_solution.pop()
            for idx in lifter_indices:
                used[idx] = False
    return False

def assign_lifters_to_boxes(boxes, lifters):
    boxes = sorted(boxes, reverse=True)
    lifters = sorted(lifters, reverse=True)
    used = [False] * len(lifters)
    best_solution = [[]]
    backtrack(boxes, lifters, 0, 8, used, [], best_solution)
    return best_solution[0]

boxes = [181, 52, 82, 130, 92, 195, 101, 106, 189, 42, 196, 46]
lifters = [45, 30, 52, 104, 30]

solution = assign_lifters_to_boxes(boxes, lifters)

output = ""
for i, step in enumerate(solution):
    output += f"Step {i+1}: {step}\n"

print(f"<<<{output}>>>")