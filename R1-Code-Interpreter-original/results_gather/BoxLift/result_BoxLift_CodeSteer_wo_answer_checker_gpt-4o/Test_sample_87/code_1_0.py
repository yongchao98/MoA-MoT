boxes = [169, 198, 165, 154, 71, 159, 55, 205, 299, 170, 122, 160, 43, 259, 246, 172]
lifters = [80, 51, 45, 141, 83, 152]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes[:]):
        # Try to find a single lifter for the box
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and lifter >= box:
                step.append((box, [j]))
                used_lifters[j] = True
                boxes.remove(box)
                break
        else:
            # Try to combine lifters to lift the box
            total_capacity = 0
            lifter_indices = []
            for j, lifter in enumerate(lifters):
                if not used_lifters[j]:
                    total_capacity += lifter
                    lifter_indices.append(j)
                    used_lifters[j] = True
                    if total_capacity >= box:
                        step.append((box, lifter_indices))
                        boxes.remove(box)
                        break
            else:
                # If we can't lift the box in this step, reset used lifters
                for idx in lifter_indices:
                    used_lifters[idx] = False

    steps.append(step)

# Print the steps
output = ""
for i, step in enumerate(steps):
    output += f"Step {i+1}: {step}\n"

print(f"<<<{output.strip()}>>>")