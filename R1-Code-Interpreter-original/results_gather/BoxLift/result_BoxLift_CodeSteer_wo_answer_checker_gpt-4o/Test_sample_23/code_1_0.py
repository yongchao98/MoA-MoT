from itertools import combinations

boxes = [254, 205, 272, 360, 320, 328, 267, 223, 337, 238, 134, 235, 388, 160, 244, 378, 349, 262, 104, 396]
lifters = [132, 86, 81, 199, 93, 92]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used_lifters):
    remaining_weight = box
    lifter_indices = []
    for i, lifter in enumerate(lifters):
        if not used_lifters[i] and remaining_weight > 0:
            if lifter <= remaining_weight:
                lifter_indices.append(i)
                remaining_weight -= lifter
                used_lifters[i] = True
    return remaining_weight <= 0, lifter_indices

def branch_and_bound(boxes, lifters, max_steps=13):
    steps = []
    used_boxes = [False] * len(boxes)
    
    def backtrack(step_count):
        if all(used_boxes):
            return True
        if step_count >= max_steps:
            return False
        
        used_lifters = [False] * len(lifters)
        step = []
        
        for i, box in enumerate(boxes):
            if used_boxes[i]:
                continue
            
            can_lift_box, lifter_indices = can_lift(box, lifters, used_lifters)
            if can_lift_box:
                step.append((box, lifter_indices))
                used_boxes[i] = True
        
        if step:
            steps.append(step)
            if backtrack(step_count + 1):
                return True
            steps.pop()
            for _, lifter_indices in step:
                for idx in lifter_indices:
                    used_lifters[idx] = False
                used_boxes[i] = False
        
        return False
    
    if backtrack(0):
        return steps
    else:
        return None

steps = branch_and_bound(boxes, lifters)

# Print the steps
if steps:
    output = "<<<"
    for idx, step in enumerate(steps):
        output += f"Step {idx + 1}: {step}\n"
    output += ">>>"
    print(output)
else:
    print("No solution found within 13 steps.")