boxes = [70, 51, 67, 51, 34, 82, 20, 93]
lifters = [64, 51, 79, 80]

def can_lift(box, lifters):
    # Check if a single lifter can lift the box
    for i, lifter in enumerate(lifters):
        if lifter >= box:
            return [i]
    # Check if a combination of lifters can lift the box
    for i in range(len(lifters)):
        for j in range(i + 1, len(lifters)):
            if lifters[i] + lifters[j] >= box:
                return [i, j]
    return None

def assign_boxes(boxes, lifters):
    steps = []
    while boxes:
        step = []
        used_lifters = set()
        remaining_boxes = []
        for box in boxes:
            assignment = can_lift(box, [l for i, l in enumerate(lifters) if i not in used_lifters])
            if assignment is not None:
                step.append((box, assignment))
                used_lifters.update(assignment)
            else:
                remaining_boxes.append(box)
        steps.append(step)
        boxes = remaining_boxes
    return steps

steps = assign_boxes(boxes, lifters)
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")