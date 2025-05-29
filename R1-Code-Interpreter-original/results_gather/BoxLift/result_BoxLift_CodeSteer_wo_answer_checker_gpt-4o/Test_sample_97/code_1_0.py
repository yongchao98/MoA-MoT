boxes = [184, 93, 275, 216, 137, 181, 31, 79, 56, 138, 81, 205, 108, 193, 230, 252]
lifters = [47, 157, 156, 45, 151, 51]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def find_combination(box, lifters, used, current_combination, index):
    if sum(current_combination) >= box:
        return current_combination
    for i in range(index, len(lifters)):
        if not used[i]:
            used[i] = True
            current_combination.append(lifters[i])
            result = find_combination(box, lifters, used, current_combination, i + 1)
            if result:
                return result
            current_combination.pop()
            used[i] = False
    return None

def assign_boxes(boxes, lifters, max_steps):
    steps = []
    for _ in range(max_steps):
        step = []
        used = [False] * len(lifters)
        for box in boxes[:]:
            lifter_combination = find_combination(box, lifters, used, [], 0)
            if lifter_combination:
                lifter_indices = [lifters.index(lifter) for lifter in lifter_combination]
                for index in lifter_indices:
                    used[index] = True
                step.append((box, lifter_indices))
                boxes.remove(box)
        steps.append(step)
        if not boxes:
            break
    return steps

steps = assign_boxes(boxes, lifters, 6)

# Format the output
output = []
for i, step in enumerate(steps):
    output.append(f"Step {i+1}: {step}")

print("<<<" + "\n".join(output) + ">>>")