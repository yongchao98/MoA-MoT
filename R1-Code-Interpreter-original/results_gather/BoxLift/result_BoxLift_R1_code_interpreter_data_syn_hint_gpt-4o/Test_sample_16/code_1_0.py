from itertools import combinations

boxes = [166, 194, 27, 144, 198, 172, 81, 117, 154, 183, 27, 184]
lifters = [97, 90, 80, 111, 76]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []

def can_lift(box, lifters, used):
    for r in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), r):
            if not any(used[i] for i in combo):
                if sum(lifters[i] for i in combo) >= box:
                    return combo
    return None

while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:
        lifter_indices = can_lift(box, lifters, used_lifters)
        if lifter_indices is not None:
            step.append((box, list(lifter_indices)))
            for i in lifter_indices:
                used_lifters[i] = True
            boxes.remove(box)
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")