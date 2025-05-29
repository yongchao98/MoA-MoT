import itertools

def find_combinations_for_weight(weight, lifters, available_indices):
    best_combination = None
    min_lifters = float('inf')
    
    for r in range(1, len(available_indices) + 1):
        for combo in itertools.combinations(available_indices, r):
            total_strength = sum(lifters[i] for i in combo)
            if total_strength >= weight and len(combo) < min_lifters:
                min_lifters = len(combo)
                best_combination = list(combo)
                if len(combo) == 1:  # Early exit if we find single lifter solution
                    return best_combination
    return best_combination

def solve_boxes(boxes, lifters):
    # Initialize
    steps = []
    remaining_boxes = list(enumerate(boxes))  # Keep track of box indices
    
    while remaining_boxes:
        current_step = []
        available_lifters = set(range(len(lifters)))
        
        # Sort remaining boxes by weight (heaviest first)
        remaining_boxes.sort(key=lambda x: x[1], reverse=True)
        
        # Try to assign boxes in current step
        i = 0
        while i < len(remaining_boxes):
            box_idx, weight = remaining_boxes[i]
            
            # Find best combination of available lifters for current box
            lifter_combo = find_combinations_for_weight(
                weight,
                lifters,
                list(available_lifters)
            )
            
            if lifter_combo:
                # Add assignment to current step
                current_step.append((weight, lifter_combo))
                # Remove used lifters from available pool
                available_lifters -= set(lifter_combo)
                # Remove assigned box
                remaining_boxes.pop(i)
            else:
                i += 1
            
            if not available_lifters:  # No more lifters available in this step
                break
        
        if current_step:
            steps.append(current_step)
        
        if len(steps) > 7:  # Check if we exceeded maximum steps
            return None
            
    return steps

# Input data
boxes = [280, 392, 142, 400, 213, 135, 295, 142, 55, 371, 151, 50, 42, 281, 283, 46, 383, 352, 204, 111]
lifters = [143, 62, 162, 128, 141, 159]

# Solve and format output
solution = solve_boxes(boxes, lifters)

if solution is None:
    print("<<<No solution found within 7 steps>>>")
else:
    output = []
    for step_idx, step in enumerate(solution, 1):
        output.append(f"Step {step_idx}: {step}")
    print("<<<" + "\n".join(output) + ">>>")