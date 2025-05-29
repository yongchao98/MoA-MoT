from itertools import combinations

def can_lift(lifters_combo, weight):
    return sum(lifters[i] for i in lifters_combo) >= weight

def find_solution():
    boxes = [21, 96, 83, 73, 100, 73, 55, 25]
    global lifters
    lifters = [62, 79, 54]
    n_lifters = len(lifters)
    boxes_left = list(enumerate(boxes))  # (index, weight)
    steps = []
    
    while boxes_left:
        step = []
        available_lifters = set(range(n_lifters))
        
        # Sort remaining boxes by weight (descending)
        boxes_left.sort(key=lambda x: x[1], reverse=True)
        
        # Try to assign boxes to lifters
        i = 0
        while i < len(boxes_left) and available_lifters:
            box_idx, weight = boxes_left[i]
            assigned = False
            
            # Try different combinations of available lifters
            for r in range(1, len(available_lifters) + 1):
                if assigned:
                    break
                for lifter_combo in combinations(available_lifters, r):
                    if can_lift(lifter_combo, weight):
                        step.append((weight, list(lifter_combo)))
                        available_lifters -= set(lifter_combo)
                        boxes_left.pop(i)
                        assigned = True
                        break
            
            if not assigned:
                i += 1
        
        steps.append(step)
        if len(steps) > 6:  # Check if we exceed maximum steps
            return None
    
    return steps

# Find and format solution
solution = find_solution()
if solution:
    output = ""
    for i, step in enumerate(solution, 1):
        output += f"Step {i}: {step}\n"
    print(output.strip())
else:
    print("No solution found within 6 steps")