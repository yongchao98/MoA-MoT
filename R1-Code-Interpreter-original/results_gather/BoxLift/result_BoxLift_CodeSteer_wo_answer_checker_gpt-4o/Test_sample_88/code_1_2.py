boxes = [88, 313, 344, 387, 392, 370, 99, 389, 340, 391, 131, 147, 334, 251, 157, 144, 290, 167, 141, 79]
lifters = [148, 193, 54, 146, 106, 53, 88]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters, used):
    # Try to find a combination of lifters that can lift the box
    def backtrack(remaining_weight, index, current_lifters):
        if remaining_weight <= 0:
            return current_lifters
        for i in range(index, len(lifters)):
            if not used[i] and lifters[i] <= remaining_weight:
                used[i] = True
                result = backtrack(remaining_weight - lifters[i], i + 1, current_lifters + [i])
                if result is not None:
                    return result
                used[i] = False
        return None

    return backtrack(box, 0, [])

steps = []
while boxes and len(steps) < 7:
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