boxes = [85, 74, 16, 69, 55, 27, 86, 63, 83, 100]
lifters = [79, 56, 71]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
max_steps = 5

while boxes and len(steps) < max_steps:
    step = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes[:]):
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and lifter >= box:
                step.append((box, [j]))
                used_lifters[j] = True
                boxes.remove(box)
                break
        else:
            # Try combinations of lifters
            for j in range(len(lifters)):
                if used_lifters[j]:
                    continue
                for k in range(j + 1, len(lifters)):
                    if used_lifters[k]:
                        continue
                    if lifters[j] + lifters[k] >= box:
                        step.append((box, [j, k]))
                        used_lifters[j] = True
                        used_lifters[k] = True
                        boxes.remove(box)
                        break
                if used_lifters[j]:
                    break
    
    steps.append(step)

# Output the steps
output = "<<<"
for i, step in enumerate(steps):
    output += f"Step {i + 1}: {step}\n"
output += ">>>"

print(output)