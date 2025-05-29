from itertools import combinations

boxes = [94, 177, 88, 131, 87, 113, 148, 86, 141, 83, 86, 73]
lifters = [55, 92, 35, 48, 96]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def find_lifter_combination(box_weight, lifters, used):
    for r in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), r):
            if all(not used[i] for i in combo):
                if sum(lifters[i] for i in combo) >= box_weight:
                    return combo
    return None

steps = []
while boxes and len(steps) < 6:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:
        combo = find_lifter_combination(box, lifters, used_lifters)
        if combo:
            step.append((box, list(combo)))
            for i in combo:
                used_lifters[i] = True
            boxes.remove(box)
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")