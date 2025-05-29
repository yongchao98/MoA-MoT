from itertools import combinations

boxes = [49, 186, 267, 243, 352, 74, 160, 115, 138, 301, 250, 145, 294, 232, 144, 293, 287, 358, 267, 266]
lifters = [60, 92, 94, 68, 72, 79]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    available_lifters = list(enumerate(lifters))
    remaining_boxes = []

    for box in boxes:
        # Try to find a combination of lifters to lift the box
        found = False
        for r in range(1, len(available_lifters) + 1):
            for combo in combinations(available_lifters, r):
                if sum(lifter[1] for lifter in combo) >= box:
                    step.append((box, [lifter[0] for lifter in combo]))
                    # Remove used lifters from available lifters
                    for lifter in combo:
                        available_lifters.remove(lifter)
                    found = True
                    break
            if found:
                break
        if not found:
            remaining_boxes.append(box)

    # Try to lift additional boxes with remaining lifters
    for box in remaining_boxes[:]:
        for r in range(1, len(available_lifters) + 1):
            for combo in combinations(available_lifters, r):
                if sum(lifter[1] for lifter in combo) >= box:
                    step.append((box, [lifter[0] for lifter in combo]))
                    for lifter in combo:
                        available_lifters.remove(lifter)
                    remaining_boxes.remove(box)
                    break

    boxes = remaining_boxes
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")