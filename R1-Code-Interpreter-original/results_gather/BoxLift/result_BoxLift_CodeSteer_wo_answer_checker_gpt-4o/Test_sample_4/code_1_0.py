boxes = [48, 293, 96, 382, 284, 61, 122, 234, 211, 109, 202, 133, 179, 293, 182, 168, 313, 101, 43, 262, 100, 231, 159, 280]
lifters = [184, 149, 64, 190, 59, 51]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used):
    # Try to find a combination of lifters that can lift the box
    for i in range(len(lifters)):
        if not used[i] and lifters[i] >= box:
            return [i]
    # If no single lifter can lift the box, try combinations
    for i in range(len(lifters)):
        if not used[i]:
            used[i] = True
            result = can_lift(box - lifters[i], lifters, used)
            used[i] = False
            if result is not None:
                return [i] + result
    return None

steps = []
while boxes and len(steps) < 8:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:
        lifter_indices = can_lift(box, lifters, used_lifters)
        if lifter_indices is not None:
            step.append((box, lifter_indices))
            boxes.remove(box)
            for idx in lifter_indices:
                used_lifters[idx] = True
    
    steps.append(step)

# Format the output
output = []
for i, step in enumerate(steps):
    output.append(f"Step {i+1}: {step}")

print("<<<" + "\n".join(output) + ">>>")