boxes = [142, 109, 62, 158, 37, 50, 180, 128, 61, 123, 44, 51]
lifters = [93, 71, 82, 91]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used):
    total_capacity = 0
    lifter_indices = []
    for i, capacity in enumerate(lifters):
        if not used[i] and total_capacity + capacity <= box:
            total_capacity += capacity
            lifter_indices.append(i)
            used[i] = True
            if total_capacity >= box:
                return True, lifter_indices
    return False, []

def find_steps(boxes, lifters, max_steps):
    steps = []
    used = [False] * len(lifters)

    for _ in range(max_steps):
        step = []
        for box in boxes[:]:
            can_lift_box, lifter_indices = can_lift(box, lifters, used)
            if can_lift_box:
                step.append((box, lifter_indices))
                boxes.remove(box)
        if step:
            steps.append(step)
            used = [False] * len(lifters)  # Reset used lifters for the next step
        if not boxes:
            break

    return steps

steps = find_steps(boxes, lifters, 5)

# Print the steps
output = ""
for i, step in enumerate(steps):
    output += f"Step {i + 1}: {step}\n"

print(f"<<<{output.strip()}>>>")