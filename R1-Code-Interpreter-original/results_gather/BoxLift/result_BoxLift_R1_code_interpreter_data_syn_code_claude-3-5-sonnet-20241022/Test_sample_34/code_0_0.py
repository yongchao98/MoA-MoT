import itertools

def find_lifting_solution():
    boxes = [399, 277, 278, 282, 209, 309, 281, 53, 325, 239, 363, 363, 308, 244, 392, 323, 255, 79, 280, 336, 395, 318, 203, 364]
    lifters = [125, 145, 57, 73, 81, 70]
    boxes.sort(reverse=True)  # Sort boxes in descending order
    
    def can_lift(weight, lifter_combo):
        return sum(lifters[i] for i in lifter_combo) >= weight
    
    steps = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(len(lifters)))
        
        while available_lifters and remaining_boxes:
            current_box = remaining_boxes[0]
            found = False
            
            # Try different combinations of available lifters
            for size in range(1, len(available_lifters) + 1):
                if found:
                    break
                for lifter_combo in itertools.combinations(available_lifters, size):
                    if can_lift(current_box, lifter_combo):
                        step.append((current_box, list(lifter_combo)))
                        available_lifters -= set(lifter_combo)
                        remaining_boxes.pop(0)
                        found = True
                        break
            
            if not found:
                break
                
        steps.append(step)
        
    # Format output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    print(result.strip())

find_lifting_solution()