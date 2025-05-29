from itertools import combinations

def can_lift(lifters, weight):
    # Try different combinations of lifters to see if they can lift the weight
    for i in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), i):
            if sum(lifters[j] for j in combo) >= weight:
                return list(combo)
    return None

def solve_box_lifting():
    boxes = [41, 69, 80, 19, 54, 85, 64, 35]
    lifter_capacities = [45, 45, 62, 74]
    boxes_left = list(enumerate(boxes))  # Keep track of original indices
    steps = []
    
    while boxes_left:
        step = []
        available_lifters = set(range(len(lifter_capacities)))
        
        # Sort boxes by weight (descending) to handle heaviest boxes first
        boxes_left.sort(key=lambda x: x[1], reverse=True)
        
        # Try to assign boxes to lifters
        i = 0
        while i < len(boxes_left) and available_lifters:
            box_idx, weight = boxes_left[i]
            
            # Find lifters that can handle this box
            available_lifter_list = [j for j in available_lifters]
            available_capacities = [lifter_capacities[j] for j in available_lifter_list]
            lifters = can_lift(available_capacities, weight)
            
            if lifters:
                # Convert local lifter indices to global indices
                global_lifters = [available_lifter_list[j] for j in lifters]
                step.append((weight, global_lifters))
                # Remove used lifters from available pool
                for lifter in global_lifters:
                    available_lifters.remove(lifter)
                # Remove assigned box
                boxes_left.pop(i)
            else:
                i += 1
        
        steps.append(step)
        
        if len(steps) > 4:  # Check if we exceed maximum steps
            return None
    
    return steps

# Solve and format output
solution = solve_box_lifting()
if solution:
    output = []
    for i, step in enumerate(solution, 1):
        output.append(f"Step {i}: {step}")
    print("<<<" + "\n".join(output) + ">>>")
else:
    print("<<<No solution found within 4 steps>>>")