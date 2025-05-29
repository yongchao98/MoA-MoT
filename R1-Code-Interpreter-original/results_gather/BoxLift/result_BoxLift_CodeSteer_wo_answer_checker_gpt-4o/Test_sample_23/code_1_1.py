boxes = [254, 205, 272, 360, 320, 328, 267, 223, 337, 238, 134, 235, 388, 160, 244, 378, 349, 262, 104, 396]
lifters = [132, 86, 81, 199, 93, 92]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used_lifters, start, current_sum):
    if current_sum >= box:
        return True
    for i in range(start, len(lifters)):
        if not used_lifters[i] and current_sum + lifters[i] <= box:
            used_lifters[i] = True
            if can_lift(box, lifters, used_lifters, i + 1, current_sum + lifters[i]):
                return True
            used_lifters[i] = False
    return False

def find_steps(boxes, lifters):
    steps = []
    used_boxes = [False] * len(boxes)
    
    while not all(used_boxes):
        step = []
        used_lifters = [False] * len(lifters)
        
        for i, box in enumerate(boxes):
            if used_boxes[i]:
                continue
            
            if can_lift(box, lifters, used_lifters, 0, 0):
                lifter_indices = [j for j, used in enumerate(used_lifters) if used]
                step.append((box, lifter_indices))
                used_boxes[i] = True
                for j in lifter_indices:
                    used_lifters[j] = True
        
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