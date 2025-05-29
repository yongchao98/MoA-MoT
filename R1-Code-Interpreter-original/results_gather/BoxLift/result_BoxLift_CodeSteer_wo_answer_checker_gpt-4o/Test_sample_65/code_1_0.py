boxes = [70, 51, 67, 51, 34, 82, 20, 93]
lifters = [64, 51, 79, 80]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used_lifters):
    # Try to find a combination of lifters that can lift the box
    lifter_indices = []
    remaining_weight = box
    
    for j, lifter in enumerate(lifters):
        if not used_lifters[j] and lifter <= remaining_weight:
            lifter_indices.append(j)
            remaining_weight -= lifter
            if remaining_weight <= 0:
                break
    
    if remaining_weight > 0:
        return None
    else:
        for index in lifter_indices:
            used_lifters[index] = True
        return lifter_indices

def assign_boxes(boxes, lifters, max_steps):
    steps = []
    used_boxes = [False] * len(boxes)
    
    def backtrack(step):
        if step >= max_steps or all(used_boxes):
            return all(used_boxes)
        
        used_lifters = [False] * len(lifters)
        step_lifts = []
        
        for i, box in enumerate(boxes):
            if used_boxes[i]:
                continue
            
            lifter_indices = can_lift(box, lifters, used_lifters)
            
            if lifter_indices is not None:
                step_lifts.append((box, lifter_indices))
                used_boxes[i] = True
        
        steps.append(step_lifts)
        
        if backtrack(step + 1):
            return True
        
        # Backtrack
        for box, lifter_indices in step_lifts:
            used_boxes[boxes.index(box)] = False
            for index in lifter_indices:
                used_lifters[index] = False
        
        steps.pop()
        return False
    
    backtrack(0)
    return steps

steps = assign_boxes(boxes, lifters, 3)

output = "<<<"
for i, step in enumerate(steps):
    output += f"Step {i+1}: {step}\n"
output += ">>>"

print(output)