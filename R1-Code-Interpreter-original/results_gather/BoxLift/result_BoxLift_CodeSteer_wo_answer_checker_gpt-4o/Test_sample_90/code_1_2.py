boxes = [71, 55, 11, 10, 71, 91, 14, 69]
lifters = [46, 54, 47]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

best_solution = []
min_steps = float('inf')
max_steps = 5

def backtrack(steps, remaining_boxes, used_lifters):
    global best_solution, min_steps

    if len(steps) >= min_steps or len(steps) >= max_steps:
        return

    if not remaining_boxes:
        if len(steps) < min_steps:
            min_steps = len(steps)
            best_solution = steps[:]
        return

    for i, box in enumerate(remaining_boxes):
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and lifter >= box:
                new_steps = steps + [(box, [j])]
                new_used_lifters = used_lifters[:]
                new_used_lifters[j] = True
                backtrack(new_steps, remaining_boxes[:i] + remaining_boxes[i+1:], new_used_lifters)
                break
        else:
            # Try to combine lifters if a single lifter can't lift the box
            for j in range(len(lifters)):
                if used_lifters[j]:
                    continue
                for k in range(j + 1, len(lifters)):
                    if used_lifters[k]:
                        continue
                    if lifters[j] + lifters[k] >= box:
                        new_steps = steps + [(box, [j, k])]
                        new_used_lifters = used_lifters[:]
                        new_used_lifters[j] = True
                        new_used_lifters[k] = True
                        backtrack(new_steps, remaining_boxes[:i] + remaining_boxes[i+1:], new_used_lifters)
                        break
                if used_lifters[j]:
                    break

# Start the backtracking process
backtrack([], boxes, [False] * len(lifters))

# Print the best solution
if best_solution:
    output = "<<<"
    for i, step in enumerate(best_solution):
        output += f"Step {i + 1}: {step}\n"
    output += ">>>"
else:
    output = "<<<No valid solution found within the constraints>>>"

print(output)