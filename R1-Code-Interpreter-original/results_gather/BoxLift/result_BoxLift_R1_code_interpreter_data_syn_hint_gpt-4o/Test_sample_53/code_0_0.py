# Box weights and lifter capacities
boxes = [108, 206, 326, 343, 45, 374, 227, 366, 110, 136, 343, 298, 201, 168, 191, 243, 365, 375, 294, 317]
lifters = [147, 80, 144, 133, 69, 200]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Function to find the minimum number of steps
def find_min_steps(boxes, lifters):
    steps = []
    while boxes:
        step = []
        remaining_lifters = lifters[:]
        for box in boxes[:]:
            # Try to lift the box with available lifters
            used_lifters = []
            total_capacity = 0
            for i, lifter in enumerate(remaining_lifters):
                if total_capacity < box:
                    used_lifters.append(i)
                    total_capacity += lifter
                if total_capacity >= box:
                    break
            if total_capacity >= box:
                # Assign lifters to this box
                step.append((box, used_lifters))
                # Remove used lifters
                for i in sorted(used_lifters, reverse=True):
                    del remaining_lifters[i]
                # Remove the box from the list
                boxes.remove(box)
        steps.append(step)
    return steps

# Get the steps
steps = find_min_steps(boxes, lifters)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")