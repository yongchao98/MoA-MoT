boxes = [100, 22, 67, 79, 32, 31, 67, 37]
lifters = [58, 78, 68, 75]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used_lifters):
    lifter_indices = []
    remaining_weight = box
    
    for i, lifter in enumerate(lifters):
        if not used_lifters[i] and lifter >= remaining_weight:
            lifter_indices.append(i)
            return lifter_indices
    
    for i, lifter in enumerate(lifters):
        if not used_lifters[i] and lifter <= remaining_weight:
            lifter_indices.append(i)
            remaining_weight -= lifter
            if remaining_weight <= 0:
                return lifter_indices
    
    return None

def backtrack(step, used_boxes, used_lifters, current_steps):
    if all(used_boxes):
        return current_steps
    
    if step >= 3:
        return None
    
    step_lifts = []
    new_used_lifters = used_lifters[:]
    
    for i, box in enumerate(boxes):
        if used_boxes[i]:
            continue
        
        lifter_indices = can_lift(box, lifters, new_used_lifters)
        
        if lifter_indices is not None:
            for idx in lifter_indices:
                new_used_lifters[idx] = True
            step_lifts.append((box, lifter_indices))
            used_boxes[i] = True
    
    if step_lifts:
        current_steps.append(step_lifts)
        result = backtrack(step + 1, used_boxes, [False] * len(lifters), current_steps)
        if result is not None:
            return result
        current_steps.pop()
    
    for i, box in enumerate(boxes):
        if used_boxes[i]:
            continue
        used_boxes[i] = False
    
    return None

used_boxes = [False] * len(boxes)
used_lifters = [False] * len(lifters)
steps = backtrack(0, used_boxes, used_lifters, [])

output = "<<<"
for i, step in enumerate(steps):
    output += f"Step {i+1}: {step}\n"
output += ">>>"

print(output)