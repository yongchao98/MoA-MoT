boxes = [50, 70, 83, 17, 21, 27, 53, 31]
lifters = [42, 53, 74, 78]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Function to find a combination of lifters to lift a box
def find_lifters_for_box(box_weight, available_lifters):
    for i in range(len(available_lifters)):
        if available_lifters[i] >= box_weight:
            return [i], available_lifters[i]
        for j in range(i + 1, len(available_lifters)):
            if available_lifters[i] + available_lifters[j] >= box_weight:
                return [i, j], available_lifters[i] + available_lifters[j]
    return [], 0

steps = []
remaining_boxes = boxes[:]
while remaining_boxes:
    step = []
    available_lifters = lifters[:]
    for box in remaining_boxes[:]:
        lifter_indices, total_capacity = find_lifters_for_box(box, available_lifters)
        if lifter_indices:
            step.append((box, lifter_indices))
            for index in sorted(lifter_indices, reverse=True):
                del available_lifters[index]
            remaining_boxes.remove(box)
    steps.append(step)

# Ensure the solution is within 3 steps
if len(steps) > 3:
    steps = "Not possible to lift all boxes in 3 steps or less"
    
print(steps)