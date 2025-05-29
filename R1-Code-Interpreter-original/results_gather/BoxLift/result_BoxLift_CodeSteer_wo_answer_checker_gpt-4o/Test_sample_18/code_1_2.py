from itertools import combinations

boxes = [161, 400, 306, 233, 179, 305, 284, 214, 288, 132, 387, 198, 274, 40, 392, 183, 275, 259, 67, 177]
lifters = [122, 116, 51, 165, 95, 150, 179]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Memoization dictionary
memo = {}

def can_lift(box_index, used_lifters):
    if box_index == len(boxes):
        return []

    if (box_index, tuple(used_lifters)) in memo:
        return memo[(box_index, tuple(used_lifters))]

    for step in range(7):
        if step not in used_lifters:
            available_lifters = [i for i in range(len(lifters)) if i not in used_lifters.get(step, [])]
            for r in range(1, len(available_lifters) + 1):
                for combo in combinations(available_lifters, r):
                    if sum(lifters[i] for i in combo) >= boxes[box_index]:
                        new_used_lifters = used_lifters.copy()
                        new_used_lifters.setdefault(step, []).extend(combo)
                        result = can_lift(box_index + 1, new_used_lifters)
                        if result is not None:
                            memo[(box_index, tuple(used_lifters))] = [(boxes[box_index], list(combo))] + result
                            return memo[(box_index, tuple(used_lifters))]

    memo[(box_index, tuple(used_lifters))] = None
    return None

solution = can_lift(0, {})

if solution is None:
    print("Unable to lift all boxes within 7 steps.")
else:
    steps = [[] for _ in range(7)]
    for box, lifter_indices in solution:
        for step in range(7):
            if all(lifter not in [l for _, lifters in steps[step] for l in lifters] for lifter in lifter_indices):
                steps[step].append((box, lifter_indices))
                break

    output = "<<<"
    for i, step in enumerate(steps):
        if step:
            output += f"Step {i + 1}: {step}\n"
    output += ">>>"
    print(output)