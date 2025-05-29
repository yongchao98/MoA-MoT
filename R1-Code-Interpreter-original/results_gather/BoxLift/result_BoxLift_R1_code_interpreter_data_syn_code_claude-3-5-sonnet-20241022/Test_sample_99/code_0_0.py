import itertools

def find_lifter_combinations(box_weight, available_lifters, lifter_capacities):
    # Try combinations of lifters from smallest to largest group size
    for i in range(1, len(available_lifters) + 1):
        for combo in itertools.combinations(available_lifters, i):
            total_capacity = sum(lifter_capacities[j] for j in combo)
            if total_capacity >= box_weight:
                return list(combo)
    return None

def solve_box_lifting():
    boxes = [289, 375, 107, 145, 257, 48, 141, 83, 136, 368, 59, 133, 186, 266, 353, 73, 66, 210, 247, 79, 342, 318, 337, 162]
    lifter_capacities = [75, 189, 162, 181, 118, 194, 137]
    
    # Sort boxes in descending order
    boxes = sorted(enumerate(boxes), key=lambda x: x[1], reverse=True)
    
    steps = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(len(lifter_capacities)))
        
        # Try to assign as many boxes as possible in current step
        for box_idx, box_weight in remaining_boxes[:]:
            if not available_lifters:
                break
                
            lifters = find_lifter_combinations(box_weight, available_lifters, lifter_capacities)
            if lifters:
                step.append((box_weight, lifters))
                available_lifters -= set(lifters)
                remaining_boxes.remove((box_idx, box_weight))
        
        if step:
            steps.append(step)
    
    # Format output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    print(result.strip())
    print(f"\nTotal steps: {len(steps)}")

solve_box_lifting()