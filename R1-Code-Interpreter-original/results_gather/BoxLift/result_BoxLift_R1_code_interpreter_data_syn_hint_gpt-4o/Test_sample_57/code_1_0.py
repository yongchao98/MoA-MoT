from itertools import combinations

boxes = [34, 92, 289, 35, 140, 174, 244, 136, 80, 252, 136, 181, 92, 235, 61, 286]
lifters = [61, 60, 150, 96, 58, 91]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes and len(steps) < 6:
    step = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes[:]):
        # Try to find the best combination of lifters for the current box
        best_combination = None
        best_capacity = 0
        
        for r in range(1, len(lifters) + 1):
            for combo in combinations(range(len(lifters)), r):
                if all(not used_lifters[j] for j in combo):
                    capacity = sum(lifters[j] for j in combo)
                    if capacity >= box and (best_combination is None or capacity < best_capacity):
                        best_combination = combo
                        best_capacity = capacity
        
        if best_combination:
            for j in best_combination:
                used_lifters[j] = True
            step.append((box, list(best_combination)))
            boxes.remove(box)
    
    steps.append(step)

for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")