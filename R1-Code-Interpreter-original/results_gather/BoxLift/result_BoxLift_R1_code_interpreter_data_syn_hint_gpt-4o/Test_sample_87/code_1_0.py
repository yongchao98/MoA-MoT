from itertools import combinations

boxes = [169, 198, 165, 154, 71, 159, 55, 205, 299, 170, 122, 160, 43, 259, 246, 172]
lifters = [80, 51, 45, 141, 83, 152]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
used_boxes = [False] * len(boxes)

def can_lift(box, lifter_indices):
    return sum(lifters[i] for i in lifter_indices) >= box

for step in range(7):
    step_lifting = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes):
        if used_boxes[i]:
            continue
        
        # Try to find the minimal combination of lifters to lift the box
        for r in range(1, len(lifters) + 1):
            for lifter_indices in combinations(range(len(lifters)), r):
                if all(not used_lifters[j] for j in lifter_indices) and can_lift(box, lifter_indices):
                    step_lifting.append((box, list(lifter_indices)))
                    used_boxes[i] = True
                    for idx in lifter_indices:
                        used_lifters[idx] = True
                    break
            if used_boxes[i]:
                break
    
    steps.append(step_lifting)
    if all(used_boxes):
        break

# Print the steps
for step_num, step in enumerate(steps, start=1):
    print(f"Step {step_num}: {step}")