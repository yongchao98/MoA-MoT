from itertools import combinations

boxes = [110, 230, 255, 47, 133, 280, 271, 275, 155, 80, 169, 89, 299, 241, 187, 234]
lifters = [160, 110, 114, 55, 54, 60]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []

def find_lifter_combination(box_weight, lifters, used_lifters):
    # Try to find the smallest combination of lifters that can lift the box
    for r in range(1, len(lifters) + 1):
        for combo in combinations(enumerate(lifters), r):
            indices, weights = zip(*combo)
            if sum(weights) >= box_weight and all(not used_lifters[i] for i in indices):
                return list(indices)
    return []

while boxes and len(steps) < 7:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:
        lifter_indices = find_lifter_combination(box, lifters, used_lifters)
        
        if lifter_indices:
            step.append((box, lifter_indices))
            for index in lifter_indices:
                used_lifters[index] = True
            boxes.remove(box)
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")