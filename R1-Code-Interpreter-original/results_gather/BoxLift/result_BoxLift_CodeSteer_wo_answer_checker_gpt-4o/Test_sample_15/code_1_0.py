boxes = [50, 70, 83, 17, 21, 27, 53, 31]
lifters = [42, 53, 74, 78]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Function to try to lift all boxes using backtracking
def lift_boxes(boxes, lifters, step, max_steps, current_steps):
    if not boxes:
        return current_steps  # All boxes are lifted

    if step >= max_steps:
        return None  # Exceeded the maximum number of steps

    available_lifters = lifters[:]
    step_lifting = []

    for box in boxes[:]:
        lifters_for_box = []
        total_capacity = 0

        for i, capacity in enumerate(available_lifters):
            if total_capacity >= box:
                break
            if capacity > 0:
                lifters_for_box.append(i)
                total_capacity += capacity
                available_lifters[i] = 0  # Mark this lifter as used in this step

        if total_capacity >= box:
            step_lifting.append((box, lifters_for_box))
            boxes.remove(box)

    if step_lifting:
        current_steps.append(step_lifting)
        result = lift_boxes(boxes, lifters, step + 1, max_steps, current_steps)
        if result is not None:
            return result
        current_steps.pop()  # Backtrack

    return None

# Try to lift all boxes in 3 steps
solution = lift_boxes(boxes[:], lifters, 0, 3, [])

# Format the output
if solution:
    output = ""
    for i, step in enumerate(solution):
        output += f"Step {i+1}: {step}\n"
    print(f"<<<{output.strip()}>>>")
else:
    print("<<<No solution found within 3 steps>>>")