from itertools import combinations

# Box weights and lifter capacities
boxes = [75, 210, 297, 242, 276, 134, 383, 358, 146, 394, 59, 350, 94, 152, 250, 136, 116, 342, 283, 87]
lifters = [101, 62, 125, 126, 138, 54]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters):
    for r in range(1, len(lifters) + 1):
        for combo in combinations(lifters, r):
            if sum(combo) >= box:
                return combo
    return None

def backtrack(boxes, lifters, steps, current_step):
    if not boxes:
        return steps

    if current_step >= 10:
        return None

    remaining_lifters = lifters[:]
    step = []
    for box in boxes[:]:
        combo = can_lift(box, remaining_lifters)
        if combo:
            step.append((box, [lifters.index(l) for l in combo]))
            for l in combo:
                remaining_lifters.remove(l)
            boxes.remove(box)

    if step:
        steps.append(step)
        result = backtrack(boxes, lifters, steps, current_step + 1)
        if result:
            return result
        steps.pop()
        boxes.extend([b[0] for b in step])

    return None

# Start backtracking
steps = backtrack(boxes, lifters, [], 0)

# Print the steps
if steps:
    output = "\n".join([f"Step {i+1}: {step}" for i, step in enumerate(steps)])
    print(f"<<<{output}>>>")
else:
    print("No solution found within 10 steps.")