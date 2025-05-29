from itertools import combinations

boxes = [289, 375, 107, 145, 257, 48, 141, 83, 136, 368, 59, 133, 186, 266, 353, 73, 66, 210, 247, 79, 342, 318, 337, 162]
lifters = [75, 189, 162, 181, 118, 194, 137]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Sort lifters in descending order
lifters.sort(reverse=True)

steps = []

def find_lifter_combination(box, lifters, used_lifters):
    for r in range(1, len(lifters) + 1):
        for combo in combinations(enumerate(lifters), r):
            indices, selected_lifters = zip(*combo)
            if sum(selected_lifters) >= box and all(not used_lifters[i] for i in indices):
                return list(indices)
    return None

while boxes and len(steps) < 7:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:  # Iterate over a copy of the list
        lifter_indices = find_lifter_combination(box, lifters, used_lifters)
        
        if lifter_indices is not None:
            step.append((box, lifter_indices))
            boxes.remove(box)
            for i in lifter_indices:
                used_lifters[i] = True
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")