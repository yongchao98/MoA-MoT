from itertools import combinations

# Box weights and lifter capacities
boxes = [108, 206, 326, 343, 45, 374, 227, 366, 110, 136, 343, 298, 201, 168, 191, 243, 365, 375, 294, 317]
lifters = [147, 80, 144, 133, 69, 200]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Function to find the minimum number of steps
def find_min_steps(boxes, lifters):
    steps = []
    while boxes:
        step = []
        used_lifters = set()
        for box in boxes[:]:
            # Try to find a single lifter or combination of lifters to lift the box
            for r in range(1, len(lifters) + 1):
                for combo in combinations(enumerate(lifters), r):
                    indices, capacities = zip(*combo)
                    if sum(capacities) >= box and not used_lifters.intersection(indices):
                        step.append((box, list(indices)))
                        used_lifters.update(indices)
                        boxes.remove(box)
                        break
                else:
                    continue
                break
        steps.append(step)
        if len(steps) > 9:
            return "Not possible within 9 steps"
    return steps

# Find the solution
solution = find_min_steps(boxes, lifters)

# Print the solution
if solution == "Not possible within 9 steps":
    print("Not possible within 9 steps")
else:
    for i, step in enumerate(solution, 1):
        print(f"Step {i}: {step}")