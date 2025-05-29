from itertools import combinations

boxes = [100, 22, 67, 79, 32, 31, 67, 37]
lifters = [58, 78, 68, 75]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used_lifters):
    for r in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), r):
            if all(not used_lifters[i] for i in combo):
                if sum(lifters[i] for i in combo) >= box:
                    return list(combo)
    return None

def backtrack(step, used_boxes, current_steps):
    if all(used_boxes):
        return current_steps
    
    if step >= 3:
        return None
    
    step_lifts = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes):
        if used_boxes[i]:
            continue
        
        lifter_indices = can_lift(box, lifters, used_lifters)
        
        if lifter_indices is not None:
            for idx in lifter_indices:
                used_lifters[idx] = True
            step_lifts.append((box, lifter_indices))
            used_boxes[i] = True
    
    if step_lifts:
        current_steps.append(step_lifts)
        result = backtrack(step + 1, used_boxes[:], current_steps)
        if result is not None:
            return result
        current_steps.pop()
    
    for i, box in enumerate(boxes):
        if used_boxes[i]:
            continue
        used_boxes[i] = False
    
    return None

used_boxes = [False] * len(boxes)
steps = backtrack(0, used_boxes, [])

if steps is None:
    print("No valid solution found within 3 steps.")
else:
    output = "<<<"
    for i, step in enumerate(steps):
        output += f"Step {i+1}: {step}\n"
    output += ">>>"
    print(output)