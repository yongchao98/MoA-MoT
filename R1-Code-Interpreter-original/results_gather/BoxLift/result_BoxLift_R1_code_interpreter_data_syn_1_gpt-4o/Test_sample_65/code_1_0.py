boxes = [70, 51, 67, 51, 34, 82, 20, 93]
lifters = [64, 51, 79, 80]

def can_lift(box, available_lifters):
    # Check if a single lifter can lift the box
    for i, lifter in enumerate(available_lifters):
        if lifter >= box:
            return [i]
    # Check if a combination of lifters can lift the box
    for i in range(len(available_lifters)):
        for j in range(i + 1, len(available_lifters)):
            if available_lifters[i] + available_lifters[j] >= box:
                return [i, j]
    return None

def assign_boxes(boxes, lifters):
    steps = []
    boxes.sort(reverse=True)  # Sort boxes by weight in descending order
    while boxes:
        step = []
        used_lifters = set()
        remaining_boxes = []
        for box in boxes:
            available_lifters = [l for i, l in enumerate(lifters) if i not in used_lifters]
            assignment = can_lift(box, available_lifters)
            if assignment is not None:
                # Map back to original lifter indices
                original_indices = [i for i, l in enumerate(lifters) if i not in used_lifters]
                step.append((box, [original_indices[i] for i in assignment]))
                used_lifters.update([original_indices[i] for i in assignment])
            else:
                remaining_boxes.append(box)
        steps.append(step)
        boxes = remaining_boxes
    return steps

steps = assign_boxes(boxes, lifters)
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")