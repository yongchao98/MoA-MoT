from itertools import combinations

# Box weights and lifter capacities
boxes = [79, 48, 16, 95, 67, 41, 62, 22]
lifters = [52, 41, 41, 67]

# Function to check if a combination of lifters can lift a box
def can_lift(box_weight, lifter_indices):
    return sum(lifters[i] for i in lifter_indices) >= box_weight

# Function to find the best assignment of lifters to boxes
def assign_lifters_to_boxes(boxes, lifters):
    steps = []
    remaining_boxes = boxes[:]
    
    for step in range(3):  # Maximum of 3 steps
        step_assignments = []
        used_lifters = set()
        
        # Try to lift each box
        for box_weight in sorted(remaining_boxes, reverse=True):
            # Find a combination of lifters that can lift this box
            for r in range(1, len(lifters) + 1):
                for lifter_indices in combinations(range(len(lifters)), r):
                    if set(lifter_indices).isdisjoint(used_lifters) and can_lift(box_weight, lifter_indices):
                        step_assignments.append((box_weight, list(lifter_indices)))
                        used_lifters.update(lifter_indices)
                        break
                else:
                    continue
                break
        
        # Remove lifted boxes from the list
        for assignment in step_assignments:
            remaining_boxes.remove(assignment[0])
        
        steps.append(step_assignments)
        
        # If all boxes are lifted, break
        if not remaining_boxes:
            break
    
    return steps

# Get the steps for lifting all boxes
steps = assign_lifters_to_boxes(boxes, lifters)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")