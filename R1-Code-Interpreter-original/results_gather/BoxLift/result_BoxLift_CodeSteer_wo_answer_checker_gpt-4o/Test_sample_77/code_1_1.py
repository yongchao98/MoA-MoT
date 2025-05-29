boxes = [253, 162, 171, 57, 68, 142, 196, 83, 94, 287, 256, 236, 153, 79, 111, 42]
lifters = [61, 62, 122, 75, 143, 131]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []

def can_lift(box, lifters, used):
    # Try to find a combination of lifters that can lift the box
    lifter_indices = []
    total_capacity = 0
    for i, lifter in enumerate(lifters):
        if not used[i] and total_capacity + lifter <= box:
            total_capacity += lifter
            lifter_indices.append(i)
            if total_capacity >= box:
                return lifter_indices
    return None

while boxes and len(steps) < 6:
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
    output += f"Step {i+1}: {step}\n"

print(f"<<<{output.strip()}>>>")