boxes = [14, 13, 59, 47, 73, 64, 58, 75]
lifters = [74, 57, 51]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used):
    # Try to find a single lifter or a combination of lifters to lift the box
    for i in range(len(lifters)):
        if not used[i] and lifters[i] >= box:
            return [i]
    for i in range(len(lifters)):
        if used[i]:
            continue
        for j in range(i + 1, len(lifters)):
            if used[j]:
                continue
            if lifters[i] + lifters[j] >= box:
                return [i, j]
    return None

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    for box in boxes[:]:
        lifter_indices = can_lift(box, lifters, used_lifters)
        if lifter_indices is not None:
            step.append((box, lifter_indices))
            for index in lifter_indices:
                used_lifters[index] = True
            boxes.remove(box)
    steps.append(step)

# Output the steps
output = ""
for i, step in enumerate(steps):
    output += f"Step {i + 1}: {step}\n"

print(f"<<<{output.strip()}>>>")