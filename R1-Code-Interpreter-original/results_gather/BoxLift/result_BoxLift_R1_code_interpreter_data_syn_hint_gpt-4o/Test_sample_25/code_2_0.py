from itertools import combinations

boxes = [17, 69, 16, 83, 95, 91, 68, 25]
lifters = [42, 68, 50, 48]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
used_boxes = [False] * len(boxes)

def find_lifter_combination(box_weight, lifters, used_lifters):
    for r in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), r):
            if all(not used_lifters[i] for i in combo):
                if sum(lifters[i] for i in combo) >= box_weight:
                    return combo
    return None

for _ in range(4):  # We have a maximum of 4 steps
    step = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes):
        if used_boxes[i]:
            continue
        
        lifter_indices = find_lifter_combination(box, lifters, used_lifters)
        
        if lifter_indices is not None:
            step.append((box, list(lifter_indices)))
            used_boxes[i] = True
            for idx in lifter_indices:
                used_lifters[idx] = True
    
    if step:
        steps.append(step)
    else:
        break

print(steps)