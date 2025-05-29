# Box weights and lifter capacities
boxes = [280, 392, 142, 400, 213, 135, 295, 142, 55, 371, 151, 50, 42, 281, 283, 46, 383, 352, 204, 111]
lifters = [143, 62, 162, 128, 141, 159]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Function to find the minimum number of steps using backtracking
def backtrack(boxes, lifters, step, steps, used_lifters):
    if not boxes:
        steps.append(step)
        return True
    if len(steps) >= 7:
        return False

    for i, box in enumerate(boxes):
        total_capacity = 0
        lifter_indices = []
        for j, lifter in enumerate(lifters):
            if not used_lifters[j]:
                total_capacity += lifter
                lifter_indices.append(j)
                if total_capacity >= box:
                    new_used_lifters = used_lifters[:]
                    for idx in lifter_indices:
                        new_used_lifters[idx] = True
                    if backtrack(boxes[:i] + boxes[i+1:], lifters, step + [(box, lifter_indices)], steps, new_used_lifters):
                        return True
                    break
    return False

# Initialize variables
steps = []
used_lifters = [False] * len(lifters)

# Start backtracking
backtrack(boxes, lifters, [], steps, used_lifters)

# Print the steps
print(steps)