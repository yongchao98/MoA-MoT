boxes = [174, 117, 74, 120, 98, 65, 112, 61, 46, 87, 176, 110]
lifters = [119, 104, 67, 98]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used):
    # Try to find a combination of lifters that can lift the box
    def backtrack(remaining_weight, start, current_lifters):
        if remaining_weight <= 0:
            return current_lifters
        for i in range(start, len(lifters)):
            if not used[i] and lifters[i] <= remaining_weight:
                used[i] = True
                result = backtrack(remaining_weight - lifters[i], i + 1, current_lifters + [i])
                if result is not None:
                    return result
                used[i] = False
        return None

    return backtrack(box, 0, [])

steps = []
while boxes and len(steps) < 5:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:
        lifter_indices = can_lift(box, lifters, used_lifters)
        if lifter_indices is not None:
            step.append((box, lifter_indices))
            boxes.remove(box)
    
    steps.append(step)

# Print the steps
output = ""
for i, step in enumerate(steps):
    output += f"Step {i + 1}: {step}\n"

print(f"<<<{output.strip()}>>>")