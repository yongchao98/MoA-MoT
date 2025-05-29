boxes = [38, 74, 55, 11, 82, 43, 34, 20]
lifters = [76, 78, 46, 63]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
remaining_boxes = boxes[:]

for step in range(3):
    step_assignments = []
    used_lifters = [False] * len(lifters)
    
    for box in remaining_boxes[:]:
        # Try to find a single lifter for the box
        for i, lifter in enumerate(lifters):
            if not used_lifters[i] and lifter >= box:
                step_assignments.append((box, [i]))
                used_lifters[i] = True
                remaining_boxes.remove(box)
                break
        else:
            # Try to find a combination of lifters for the box
            for i in range(len(lifters)):
                if used_lifters[i]:
                    continue
                for j in range(i + 1, len(lifters)):
                    if used_lifters[j]:
                        continue
                    if lifters[i] + lifters[j] >= box:
                        step_assignments.append((box, [i, j]))
                        used_lifters[i] = True
                        used_lifters[j] = True
                        remaining_boxes.remove(box)
                        break
                if box not in remaining_boxes:
                    break
    
    steps.append(step_assignments)
    if not remaining_boxes:
        break

# Print the steps
output = ""
for i, step in enumerate(steps):
    output += f"Step {i + 1}: {step}\n"

print(f"<<<{output.strip()}>>>")