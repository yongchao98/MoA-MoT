from itertools import combinations

def can_lift(lifters, weight):
    # Check all possible combinations of lifters
    for r in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), r):
            if sum(lifters[i] for i in combo) >= weight:
                return list(combo)
    return None

def solve_box_lifting():
    boxes = [215, 113, 62, 255, 246, 148, 283, 149, 295, 250, 114, 48, 245, 145, 153, 50]
    lifters = [82, 55, 159, 56, 88]
    
    # Sort boxes in descending order
    boxes = sorted([(w, i) for i, w in enumerate(boxes)], reverse=True)
    solution = []
    
    while boxes:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # Try to assign as many boxes as possible in current step
        i = 0
        while i < len(boxes):
            weight, orig_idx = boxes[i]
            
            # Find valid lifter combination
            valid_lifters = [lifters[j] for j in available_lifters]
            lifter_indices = can_lift([lifters[j] for j in available_lifters], weight)
            
            if lifter_indices is not None:
                # Convert relative indices to absolute indices
                abs_indices = [list(available_lifters)[j] for j in lifter_indices]
                step.append((weight, abs_indices))
                # Remove used lifters
                available_lifters -= set(abs_indices)
                boxes.pop(i)
            else:
                i += 1
                
        if step:
            solution.append(step)
        
    # Print solution
    for i, step in enumerate(solution, 1):
        print(f"Step {i}: {step}")
    print(f"\nTotal steps: {len(solution)}")

solve_box_lifting()