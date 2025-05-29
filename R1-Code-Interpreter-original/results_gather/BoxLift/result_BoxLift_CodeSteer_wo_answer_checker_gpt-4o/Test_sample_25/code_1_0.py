from itertools import combinations

def can_lift(box, lifters, used):
    for r in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), r):
            if all(not used[i] for i in combo) and sum(lifters[i] for i in combo) >= box:
                for i in combo:
                    used[i] = True
                return list(combo)
    return None

def branch_and_bound(boxes, lifters, step, steps, used):
    if not boxes:
        return True
    
    if step >= 4:
        return False
    
    current_step = []
    remaining_boxes = boxes[:]
    
    for box in boxes:
        lifter_indices = can_lift(box, lifters, used)
        if lifter_indices is not None:
            current_step.append((box, lifter_indices))
            remaining_boxes.remove(box)
    
    if current_step:
        steps.append(current_step)
        if branch_and_bound(remaining_boxes, lifters, step + 1, steps, [False] * len(lifters)):
            return True
        steps.pop()
    
    return False

def assign_lifters_to_boxes(boxes, lifters):
    boxes.sort(reverse=True)
    lifters.sort(reverse=True)
    steps = []
    used = [False] * len(lifters)
    
    if branch_and_bound(boxes, lifters, 0, steps, used):
        return steps
    else:
        return "Not possible to lift all boxes in 4 steps or less."

boxes = [17, 69, 16, 83, 95, 91, 68, 25]
lifters = [42, 68, 50, 48]

steps = assign_lifters_to_boxes(boxes, lifters)

output = "<<<"
for i, step in enumerate(steps):
    output += f"Step {i+1}: {step}\n"
output += ">>>"

print(output)