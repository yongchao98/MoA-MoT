boxes = [202, 297, 200, 241, 156, 177, 70, 78, 133, 176, 75, 70, 167, 240, 299, 91]
lifters = [56, 160, 118, 138, 160, 150]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def find_combination(box, lifters, used):
    # Try to find a combination of lifters to lift the box
    combination = []
    total_capacity = 0
    for i, lifter in enumerate(lifters):
        if not used[i]:
            combination.append(i)
            total_capacity += lifter
            if total_capacity >= box:
                return combination
    return None

steps = []
while boxes and len(steps) < 5:
    step = []
    used_lifters = [False] * len(lifters)
    for box in boxes[:]:
        combination = find_combination(box, lifters, used_lifters)
        if combination:
            step.append((box, combination))
            for index in combination:
                used_lifters[index] = True
            boxes.remove(box)
    steps.append(step)

# Format the output
output = ""
for i, step in enumerate(steps):
    output += f"Step {i + 1}: {step}\n"

print(f"<<<{output.strip()}>>>")