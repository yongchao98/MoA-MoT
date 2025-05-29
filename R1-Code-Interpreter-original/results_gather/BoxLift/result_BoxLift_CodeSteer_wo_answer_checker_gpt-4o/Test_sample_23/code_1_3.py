from heapq import heappush, heappop

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

def find_steps(boxes, lifters):
    steps = []
    used_boxes = [False] * len(boxes)
    
    while not all(used_boxes):
        step = []
        used_lifters = [False] * len(lifters)
        
        for i, box in enumerate(boxes):
            if used_boxes[i]:
                continue
            
            can_lift_box, lifter_indices = can_lift(box, lifters, used_lifters)
            if can_lift_box:
                step.append((box, lifter_indices))
                used_boxes[i] = True
        
        steps.append(step)
        if len(steps) > 13:
            return None  # Exceeded step limit
    
    return steps

steps = find_steps(boxes, lifters)

# Print the steps
if steps:
    output = "<<<"
    for idx, step in enumerate(steps):
        output += f"Step {idx + 1}: {step}\n"
    output += ">>>"
    print(output)
else:
    print("No solution found within 13 steps.")