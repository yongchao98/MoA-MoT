from itertools import combinations

def can_lift(lifters, weight):
    # Check all possible combinations of lifters
    for i in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), i):
            if sum(lifters[j] for j in combo) >= weight:
                return list(combo)
    return None

def solve_box_lifting():
    boxes = [38, 74, 55, 11, 82, 43, 34, 20]
    lifter_capacities = [76, 78, 46, 63]
    solution = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes and len(solution) < 3:
        step = []
        available_lifters = set(range(len(lifter_capacities)))
        
        # Sort boxes in descending order
        boxes_with_indices = [(w, i) for i, w in enumerate(remaining_boxes)]
        boxes_with_indices.sort(reverse=True)
        
        new_remaining = remaining_boxes.copy()
        
        for weight, orig_idx in boxes_with_indices:
            if not available_lifters:
                break
                
            # Convert available_lifters to list for combination checking
            curr_lifters = [lifter_capacities[i] for i in available_lifters]
            lifter_indices = list(available_lifters)
            lifters_needed = can_lift([curr_lifters[i] for i in range(len(curr_lifters))], weight)
            
            if lifters_needed:
                # Convert to actual lifter indices
                actual_lifters = [lifter_indices[i] for i in lifters_needed]
                step.append((weight, actual_lifters))
                new_remaining.remove(weight)
                
                # Remove used lifters
                for lifter in actual_lifters:
                    available_lifters.remove(lifter)
        
        if step:
            solution.append(step)
            remaining_boxes = new_remaining
        else:
            break
    
    # Format and print solution
    if not remaining_boxes and len(solution) <= 3:
        result = ""
        for i, step in enumerate(solution, 1):
            result += f"Step {i}: {step}\n"
        print(result.strip())
    else:
        print("No solution found within 3 steps")

solve_box_lifting()