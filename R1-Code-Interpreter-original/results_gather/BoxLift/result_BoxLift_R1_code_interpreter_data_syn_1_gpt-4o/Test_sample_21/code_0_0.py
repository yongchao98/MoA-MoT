# Define the weights of the boxes and the capacities of the lifters
boxes = [38, 97, 21, 92, 195, 187, 95, 140, 37, 85, 25, 165]
lifters = [96, 32, 88, 48]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Sort lifters in descending order with their indices
lifters_with_indices = sorted([(capacity, index) for index, capacity in enumerate(lifters)], reverse=True)

# Initialize the steps
steps = []

# Function to find a combination of lifters that can lift a given box
def find_lifters_for_box(box_weight, available_lifters):
    used_lifters = []
    total_capacity = 0
    for capacity, index in available_lifters:
        if total_capacity < box_weight:
            used_lifters.append(index)
            total_capacity += capacity
        if total_capacity >= box_weight:
            break
    if total_capacity >= box_weight:
        return used_lifters
    return None

# Try to lift all boxes in 6 steps or less
for step in range(6):
    step_lifting = []
    available_lifters = lifters_with_indices[:]
    remaining_boxes = []
    
    for box in boxes:
        lifters_for_box = find_lifters_for_box(box, available_lifters)
        if lifters_for_box is not None:
            step_lifting.append((box, lifters_for_box))
            # Remove used lifters from available lifters
            available_lifters = [(cap, idx) for cap, idx in available_lifters if idx not in lifters_for_box]
        else:
            remaining_boxes.append(box)
    
    steps.append(step_lifting)
    boxes = remaining_boxes
    
    if not boxes:
        break

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")