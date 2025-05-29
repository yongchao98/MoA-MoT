from itertools import combinations

def can_lift(lifters, weight):
    # Check all possible combinations of lifters to see if they can lift the weight
    for i in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), i):
            if sum(lifters[j] for j in combo) >= weight:
                return list(combo)
    return None

def solve_box_lifting():
    boxes = [76, 97, 14, 30, 36, 83, 45, 46]
    lifter_capacities = [69, 59, 68, 50]
    boxes_with_index = list(enumerate(boxes))
    solution = []
    
    while boxes_with_index:
        step = []
        available_lifters = set(range(len(lifter_capacities)))
        
        # Sort boxes by weight (descending) for better optimization
        boxes_with_index.sort(key=lambda x: x[1], reverse=True)
        
        # Try to assign boxes to lifters
        i = 0
        while i < len(boxes_with_index) and available_lifters:
            box_idx, weight = boxes_with_index[i]
            
            # Get available lifter capacities
            current_lifters = [lifter_capacities[j] for j in available_lifters]
            lifter_indices = [j for j in available_lifters]
            
            # Try to find lifters for this box
            required_lifters = can_lift([current_lifters[i] for i in range(len(current_lifters))], weight)
            
            if required_lifters is not None:
                # Convert local lifter indices to global indices
                global_lifters = [lifter_indices[j] for j in required_lifters]
                step.append((weight, global_lifters))
                # Remove used lifters from available pool
                for lifter in global_lifters:
                    available_lifters.remove(lifter)
                # Remove assigned box
                boxes_with_index.pop(i)
            else:
                i += 1
        
        if step:
            solution.append(step)
    
    # Format the solution
    result = ""
    for i, step in enumerate(solution, 1):
        result += f"Step {i}: {step}\n"
    
    print(result.strip())

solve_box_lifting()