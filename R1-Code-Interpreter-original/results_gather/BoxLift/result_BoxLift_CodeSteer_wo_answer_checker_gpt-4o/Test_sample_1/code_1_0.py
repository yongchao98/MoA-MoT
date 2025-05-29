boxes = [71, 58, 105, 85, 123, 125, 84, 105, 141, 73, 102, 46]
lifters = [77, 45, 120, 109]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
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
            # Try to combine lifters if a single lifter can't lift the box
            for j in range(len(lifters)):
                if used_lifters[j]:
                    continue
                total_capacity = lifters[j]
                lifter_indices = [j]
                for k in range(j + 1, len(lifters)):
                    if used_lifters[k]:
                        continue
                    if total_capacity >= box:
                        break
                    total_capacity += lifters[k]
                    lifter_indices.append(k)
                if total_capacity >= box:
                    for idx in lifter_indices:
                        used_lifters[idx] = True
                    step.append((box, lifter_indices))
                    boxes.remove(box)
                    break
    
    steps.append(step)

# Print the steps
output = ""
for i, step in enumerate(steps):
    output += f"Step {i + 1}: {step}\n"

print(f"<<<{output.strip()}>>>")