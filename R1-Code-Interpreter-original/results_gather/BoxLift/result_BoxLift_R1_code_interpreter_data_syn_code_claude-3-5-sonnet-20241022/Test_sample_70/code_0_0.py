from itertools import combinations

def can_lift(lifters_combo, weight):
    return sum(lifters_combo) >= weight

def find_solution():
    boxes = [82, 61, 91, 29, 20, 43, 57, 79]
    lifters = [49, 80, 79]
    boxes_left = boxes.copy()
    steps = []
    
    while boxes_left and len(steps) < 4:
        step = []
        available_lifters = set(range(len(lifters)))
        boxes_left_sorted = sorted(boxes_left, reverse=True)
        
        for box in boxes_left_sorted:
            if not available_lifters:
                break
                
            # Try all possible combinations of available lifters
            for r in range(1, len(available_lifters) + 1):
                found = False
                for lifter_indices in combinations(available_lifters, r):
                    lifter_capacities = [lifters[i] for i in lifter_indices]
                    if sum(lifter_capacities) >= box:
                        step.append((box, list(lifter_indices)))
                        available_lifters -= set(lifter_indices)
                        boxes_left.remove(box)
                        found = True
                        break
                if found:
                    break
        
        if step:
            steps.append(step)
        else:
            return None
    
    return steps if not boxes_left else None

# Find and format solution
solution = find_solution()
if solution:
    result = ""
    for i, step in enumerate(solution, 1):
        result += f"Step {i}: {step}\n"
    print(result.strip())
else:
    print("No solution found within 4 steps")