from itertools import combinations

boxes = [159, 180, 64, 184, 110, 43, 128, 47, 47, 138, 62, 97]
lifters = [81, 98, 62, 46]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Function to find a combination of lifters that can lift a given box
def find_lifters_for_box(box, lifters):
    for r in range(1, len(lifters) + 1):
        for combo in combinations(enumerate(lifters), r):
            indices, weights = zip(*combo)
            if sum(weights) >= box:
                return list(indices)
    return []

steps = []
while boxes and len(steps) < 7:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:  # Iterate over a copy of the list
        lifter_indices = find_lifters_for_box(box, [l for i, l in enumerate(lifters) if not used_lifters[i]])
        
        if lifter_indices:
            step.append((box, lifter_indices))
            for index in lifter_indices:
                used_lifters[index] = True
            boxes.remove(box)
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")