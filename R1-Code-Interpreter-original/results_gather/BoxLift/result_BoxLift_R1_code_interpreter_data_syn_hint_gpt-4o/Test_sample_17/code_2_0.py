from itertools import combinations

boxes = [130, 182, 34, 60, 202, 148, 160, 265, 216, 93, 133, 245, 165, 104, 143, 220]
lifters = [76, 78, 96, 122, 74, 80]

# Sort boxes in descending order
boxes.sort(reverse=True)

def can_lift(box, lifters):
    for r in range(1, len(lifters) + 1):
        for combo in combinations(lifters, r):
            if sum(combo) >= box:
                return list(combo)
    return None

def lift_boxes(boxes, lifters, steps, max_steps):
    if not boxes:
        return steps
    if len(steps) >= max_steps:
        return None
    
    for i, box in enumerate(boxes):
        lifter_combo = can_lift(box, lifters)
        if lifter_combo:
            new_lifters = [l for l in lifters if l not in lifter_combo]
            new_steps = steps + [(box, lifter_combo)]
            result = lift_boxes(boxes[:i] + boxes[i+1:], new_lifters, new_steps, max_steps)
            if result is not None:
                return result
    return None

result = lift_boxes(boxes, lifters, [], 7)

# Print the result
if result:
    for i, step in enumerate(result):
        print(f"Step {i+1}: {step}")
else:
    print("Could not lift all boxes within 7 steps.")