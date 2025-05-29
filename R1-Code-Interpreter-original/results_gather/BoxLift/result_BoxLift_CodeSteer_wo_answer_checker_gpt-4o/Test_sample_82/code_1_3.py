from itertools import combinations

boxes = [99, 159, 55, 86, 215, 56, 271, 167, 212, 87, 78, 211, 246, 153, 227, 287]
lifters = [143, 50, 70, 79, 136]

# Sort boxes in descending order
boxes.sort(reverse=True)

def find_lifter_combinations(lifters, box_weight):
    # Find all combinations of lifters that can lift the box
    for r in range(1, len(lifters) + 1):
        for combo in combinations(enumerate(lifters), r):
            indices, weights = zip(*combo)
            if sum(weights) >= box_weight:
                return list(indices)
    return None

def solve_box_lift(boxes, lifters, max_steps):
    steps = []
    while boxes and len(steps) < max_steps:
        step = []
        used = [False] * len(lifters)
        for box in boxes[:]:
            lifter_indices = find_lifter_combinations([l for i, l in enumerate(lifters) if not used[i]], box)
            if lifter_indices is not None:
                step.append((box, lifter_indices))
                boxes.remove(box)
                for index in lifter_indices:
                    used[index] = True
        if step:
            steps.append(step)
        else:
            break

    if not boxes:
        output = []
        for i, step in enumerate(steps):
            output.append(f"Step {i + 1}: {step}")
        print("<<<" + "\n".join(output) + ">>>")
    else:
        print("No solution found within the given steps.")

solve_box_lift(boxes, lifters, 7)