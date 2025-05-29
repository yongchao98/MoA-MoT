from itertools import combinations

boxes = [159, 197, 211, 302, 51, 251, 299, 229, 384, 196, 377, 283, 127, 312, 228, 280, 153, 199, 169, 151]
lifters = [52, 80, 175, 75, 136, 143]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes and len(steps) < 10:
    step = []
    used_lifters = [False] * len(lifters)
    remaining_boxes = []
    
    for box in boxes:
        lifter_indices = []
        for r in range(1, len(lifters) + 1):
            for combo in combinations(enumerate(lifters), r):
                indices, capacities = zip(*combo)
                if sum(capacities) >= box and all(not used_lifters[i] for i in indices):
                    lifter_indices = list(indices)
                    for i in indices:
                        used_lifters[i] = True
                    break
            if lifter_indices:
                break
        
        if lifter_indices:
            step.append((box, lifter_indices))
        else:
            remaining_boxes.append(box)
    
    steps.append(step)
    boxes = remaining_boxes

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")