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

def assign_boxes_to_steps(boxes, lifters, max_steps):
    steps = []
    remaining_boxes = boxes[:]
    while remaining_boxes and len(steps) < max_steps:
        step = []
        available_lifters = lifters[:]
        used_lifters = set()
        for box in remaining_boxes[:]:
            lifter_indices, total_capacity = find_lifters_for_box(box, available_lifters)
            if lifter_indices and not any(index in used_lifters for index in lifter_indices):
                step.append((box, lifter_indices))
                used_lifters.update(lifter_indices)
                for index in sorted(lifter_indices, reverse=True):
                    del available_lifters[index]
                remaining_boxes.remove(box)
        steps.append(step)
    return steps if not remaining_boxes else "Not possible to lift all boxes in 3 steps or less"

steps = assign_boxes_to_steps(boxes, lifters, 3)
print(steps)