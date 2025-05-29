from itertools import combinations

boxes = [75, 298, 70, 200, 278, 35, 278, 56, 270, 212, 161, 98, 244, 292, 177, 145]
lifters = [105, 102, 125, 93, 123]

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
                if combo not in used_lifters:
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