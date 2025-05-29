from itertools import combinations

boxes = [281, 159, 239, 61, 219, 44, 225, 192, 74, 369, 318, 366, 62, 186, 373, 243, 131, 172, 66, 167, 400, 330, 351, 245]
lifters = [164, 182, 169, 127, 118, 60]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes:
    step = []
    remaining_lifters = lifters[:]
    lifter_indices = list(range(len(lifters)))
    used_lifters = set()
    
    for box in boxes[:]:
        # Try to find a combination of lifters to lift the box
        for r in range(1, len(lifters) + 1):
            for combo in combinations(lifter_indices, r):
                if not any(i in used_lifters for i in combo):
                    total_capacity = sum(remaining_lifters[i] for i in combo)
                    if total_capacity >= box:
                        step.append((box, list(combo)))
                        used_lifters.update(combo)
                        boxes.remove(box)
                        break
            if box not in boxes:
                break
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")