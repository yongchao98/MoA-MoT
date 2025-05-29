# Define the weights of the boxes and the capacities of the lifters
boxes = [17, 69, 16, 83, 95, 91, 68, 25]
lifters = [42, 68, 50, 48]

# Sort boxes in descending order and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Function to find the best combination of lifters to lift a box
def find_lifters_for_box(box_weight, available_lifters):
    # Try to find a single lifter first
    for i, lifter in enumerate(available_lifters):
        if lifter >= box_weight:
            return [i]
    
    # If no single lifter can lift the box, try combinations
    from itertools import combinations
    for r in range(2, len(available_lifters) + 1):
        for combo in combinations(enumerate(available_lifters), r):
            indices, lifter_weights = zip(*combo)
            if sum(lifter_weights) >= box_weight:
                return list(indices)
    
    return None

# Assign boxes to lifters in steps
steps = []
while boxes and len(steps) < 4:
    step = []
    available_lifters = lifters[:]
    used_lifters = set()
    
    for box in boxes[:]:
        lifter_indices = find_lifters_for_box(box, available_lifters)
        if lifter_indices is not None:
            step.append((box, lifter_indices))
            boxes.remove(box)
            for index in lifter_indices:
                used_lifters.add(index)
            # Update available lifters for this step
            available_lifters = [lifter for i, lifter in enumerate(lifters) if i not in used_lifters]
    
    steps.append(step)

# Print the steps
print(steps)