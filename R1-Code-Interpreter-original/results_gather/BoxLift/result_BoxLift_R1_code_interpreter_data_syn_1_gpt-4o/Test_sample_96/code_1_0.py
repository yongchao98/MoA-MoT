from itertools import combinations

boxes = [169, 106, 153, 20, 138, 80, 94, 190, 142, 196, 186, 87]
lifters = [55, 85, 43, 100, 47]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []

def find_lifter_combination(box_weight, lifters):
    # Try all combinations of lifters to find a valid one
    for r in range(1, len(lifters) + 1):
        for combo in combinations(enumerate(lifters), r):
            indices, weights = zip(*combo)
            if sum(weights) >= box_weight:
                return list(indices)
    return None

while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:
        lifter_indices = find_lifter_combination(box, [l for i, l in enumerate(lifters) if not used_lifters[i]])
        
        if lifter_indices is not None:
            step.append((box, lifter_indices))
            boxes.remove(box)
            for index in lifter_indices:
                used_lifters[index] = True
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")