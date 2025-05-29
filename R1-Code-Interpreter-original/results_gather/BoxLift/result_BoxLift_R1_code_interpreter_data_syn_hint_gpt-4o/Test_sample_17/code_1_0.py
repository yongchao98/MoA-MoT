from itertools import combinations

boxes = [130, 182, 34, 60, 202, 148, 160, 265, 216, 93, 133, 245, 165, 104, 143, 220]
lifters = [76, 78, 96, 122, 74, 80]

# Sort boxes in descending order
boxes.sort(reverse=True)

def find_lifter_combination(box, lifters):
    for r in range(1, len(lifters) + 1):
        for combo in combinations(lifters, r):
            if sum(combo) >= box:
                return list(combo)
    return None

steps = []
while boxes and len(steps) < 7:
    step = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes):
        lifter_combo = find_lifter_combination(box, [l for j, l in enumerate(lifters) if not used_lifters[j]])
        if lifter_combo:
            step.append((box, lifter_combo))
            for l in lifter_combo:
                used_lifters[lifters.index(l)] = True
    
    if step:
        steps.append(step)
        for box, _ in step:
            boxes.remove(box)
    else:
        break

# Print the steps
if boxes:
    print("Could not lift all boxes within 7 steps.")
else:
    for i, step in enumerate(steps):
        print(f"Step {i+1}: {step}")