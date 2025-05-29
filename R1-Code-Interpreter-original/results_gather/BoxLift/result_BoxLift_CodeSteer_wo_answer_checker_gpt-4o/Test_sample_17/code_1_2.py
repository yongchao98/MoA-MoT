boxes = [130, 182, 34, 60, 202, 148, 160, 265, 216, 93, 133, 245, 165, 104, 143, 220]
lifters = [76, 78, 96, 122, 74, 80]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def branch_and_bound(boxes, lifters, current_step, steps, best_steps):
    if not boxes:
        if len(steps) < len(best_steps[0]):
            best_steps[0] = steps[:]
        return
    
    if len(steps) >= len(best_steps[0]):
        return
    
    used_lifters = [False] * len(lifters)
    step = []
    
    for box in boxes[:]:
        total_capacity = 0
        lifter_indices = []
        
        for i, lifter in enumerate(lifters):
            if not used_lifters[i] and total_capacity + lifter <= box:
                total_capacity += lifter
                lifter_indices.append(i)
                used_lifters[i] = True
                if total_capacity >= box:
                    break
        
        if total_capacity >= box:
            step.append((box, lifter_indices))
            boxes.remove(box)
    
    if step:
        steps.append(step)
        branch_and_bound(boxes, lifters, current_step + 1, steps, best_steps)
        steps.pop()
        boxes.extend([b[0] for b in step])

best_steps = [[(None, [])] * 8]  # Start with a large number of steps
branch_and_bound(boxes, lifters, 0, [], best_steps)

# Print the best steps
output = "<<<"
for i, step in enumerate(best_steps[0]):
    output += f"Step {i+1}: {step}\n"
output += ">>>"

print(output)