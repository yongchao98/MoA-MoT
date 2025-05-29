# Define the weights of the boxes and the capacities of the lifters
boxes = [57, 112, 121, 49, 166, 143, 169, 43, 90, 66, 23, 143]
lifters = [103, 49, 95, 97, 98]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Function to find the best combination of lifters for a given box
def find_lifters_for_box(box_weight, available_lifters):
    from itertools import combinations
    for r in range(1, len(available_lifters) + 1):
        for combo in combinations(available_lifters, r):
            if sum(combo) >= box_weight:
                return list(combo)
    return []

# Initialize steps
steps = []

# Try to lift all boxes in 4 steps
for step in range(4):
    step_lifting = []
    used_lifters = set()
    
    for box in boxes[:]:
        available_lifters = [l for i, l in enumerate(lifters) if i not in used_lifters]
        lifter_combo = find_lifters_for_box(box, available_lifters)
        
        if lifter_combo:
            step_lifting.append((box, [lifters.index(l) for l in lifter_combo]))
            boxes.remove(box)
            for l in lifter_combo:
                used_lifters.add(lifters.index(l))
    
    steps.append(step_lifting)
    
    if not boxes:
        break

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")