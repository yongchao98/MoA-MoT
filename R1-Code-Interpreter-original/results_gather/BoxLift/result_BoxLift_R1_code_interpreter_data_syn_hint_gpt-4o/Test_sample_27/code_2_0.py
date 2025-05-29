from itertools import combinations

# Define the weights of the boxes and the capacities of the lifters
boxes = [110, 230, 255, 47, 133, 280, 271, 275, 155, 80, 169, 89, 299, 241, 187, 234]
lifters = [160, 110, 114, 55, 54, 60]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Initialize the steps
steps = []

# Function to find a combination of lifters to lift a box
def find_lifters_for_box(box_weight, available_lifters):
    for r in range(1, len(available_lifters) + 1):
        for indices in combinations(range(len(available_lifters)), r):
            total_capacity = sum(available_lifters[i] for i in indices)
            if total_capacity >= box_weight:
                return indices
    return None

# Try to lift all boxes in 7 steps or less
for step in range(7):
    step_lifting = []
    available_lifters = lifters[:]
    used_lifters = set()
    
    for box in boxes[:]:
        lifter_indices = find_lifters_for_box(box, available_lifters)
        if lifter_indices is not None:
            step_lifting.append((box, list(lifter_indices)))
            boxes.remove(box)
            for index in lifter_indices:
                used_lifters.add(index)
    
    # Remove used lifters from available lifters
    available_lifters = [capacity for i, capacity in enumerate(lifters) if i not in used_lifters]
    
    steps.append(step_lifting)
    
    if not boxes:
        break

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")