boxes = [253, 162, 171, 57, 68, 142, 196, 83, 94, 287, 256, 236, 153, 79, 111, 42]
lifters = [61, 62, 122, 75, 143, 131]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []

def find_combination(box, lifters, used, start=0, current_combination=[], current_sum=0):
    if current_sum >= box:
        return current_combination
    for i in range(start, len(lifters)):
        if not used[i] and current_sum + lifters[i] <= box:
            used[i] = True
            result = find_combination(box, lifters, used, i + 1, current_combination + [i], current_sum + lifters[i])
            if result is not None:
                return result
            used[i] = False
    return None

while boxes and len(steps) < 6:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:
        lifter_indices = find_combination(box, lifters, used_lifters)
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