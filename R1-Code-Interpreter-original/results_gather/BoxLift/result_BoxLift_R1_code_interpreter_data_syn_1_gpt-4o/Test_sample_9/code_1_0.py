from itertools import combinations

boxes = [355, 40, 260, 256, 250, 107, 301, 298, 72, 316, 50, 319, 243, 51, 393, 268, 141, 178, 270, 158, 354, 336, 55, 84]
lifters = [100, 90, 66, 141, 184, 78, 50]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:
        # Try to find a combination of lifters to lift the box
        lifter_indices = []
        remaining_weight = box
        
        # Check all combinations of lifters
        for r in range(1, len(lifters) + 1):
            for combo in combinations(enumerate(lifters), r):
                indices, weights = zip(*combo)
                if not any(used_lifters[i] for i in indices) and sum(weights) >= remaining_weight:
                    lifter_indices = list(indices)
                    for i in indices:
                        used_lifters[i] = True
                    break
            if lifter_indices:
                break
        
        if lifter_indices:
            step.append((box, lifter_indices))
            boxes.remove(box)
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")