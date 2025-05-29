import itertools

def find_lifting_solution():
    boxes = [399, 277, 278, 282, 209, 309, 281, 53, 325, 239, 363, 363, 308, 244, 392, 323, 255, 79, 280, 336, 395, 318, 203, 364]
    lifters = [125, 145, 57, 73, 81, 70]
    boxes.sort(reverse=True)
    
    def can_lift(weight, lifter_combo):
        return sum(lifters[i] for i in lifter_combo) >= weight
    
    def find_best_lifter_combo(weight, available_lifters):
        for size in range(1, len(available_lifters) + 1):
            for combo in itertools.combinations(available_lifters, size):
                if can_lift(weight, combo):
                    return list(combo)
        return None
    
    steps = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # Try to handle heavy boxes first
        i = 0
        while i < len(remaining_boxes) and available_lifters:
            current_box = remaining_boxes[i]
            lifter_combo = find_best_lifter_combo(current_box, available_lifters)
            
            if lifter_combo:
                step.append((current_box, lifter_combo))
                available_lifters -= set(lifter_combo)
                remaining_boxes.pop(i)
            else:
                i += 1
            
            # Try to fit small boxes with remaining lifters
            if available_lifters:
                for j in range(len(remaining_boxes)-1, -1, -1):
                    if j != i-1:
                        small_box = remaining_boxes[j]
                        if small_box < sum(lifters[k] for k in available_lifters):
                            lifter_combo = find_best_lifter_combo(small_box, available_lifters)
                            if lifter_combo:
                                step.append((small_box, lifter_combo))
                                available_lifters -= set(lifter_combo)
                                remaining_boxes.pop(j)
                                if j < i:
                                    i -= 1
        
        if step:
            steps.append(step)
        
    # Format output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    print(result.strip())

find_lifting_solution()