from itertools import combinations

def can_lift(lifters, weight):
    # Try all possible combinations of lifters
    for r in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), r):
            if sum(lifters[i] for i in combo) >= weight:
                return list(combo)
    return None

def solve_box_lifting():
    boxes = [169, 106, 153, 20, 138, 80, 94, 190, 142, 196, 186, 87]
    lifter_capacities = [55, 85, 43, 100, 47]
    
    # Sort boxes in descending order
    boxes = sorted(boxes, reverse=True)
    
    solution = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(len(lifter_capacities)))
        
        # Try to assign boxes to lifters
        for box in remaining_boxes[:]:
            if not available_lifters:
                break
                
            # Find lifters that can lift this box
            lifters_subset = [lifter_capacities[i] for i in available_lifters]
            lifter_indices = can_lift([lifter_capacities[i] for i in available_lifters], box)
            
            if lifter_indices is not None:
                # Convert relative indices to absolute indices
                abs_indices = [list(available_lifters)[i] for i in lifter_indices]
                step.append((box, abs_indices))
                remaining_boxes.remove(box)
                # Remove used lifters
                for idx in abs_indices:
                    available_lifters.remove(idx)
        
        if step:
            solution.append(step)
        else:
            break
    
    # Print solution
    result = ""
    for i, step in enumerate(solution, 1):
        result += f"Step {i}: {step}\n"
    
    print(result.strip())

solve_box_lifting()