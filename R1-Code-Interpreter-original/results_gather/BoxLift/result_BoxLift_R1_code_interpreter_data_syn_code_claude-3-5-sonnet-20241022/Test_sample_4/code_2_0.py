import itertools

def find_lifting_solution():
    boxes = [48, 293, 96, 382, 284, 61, 122, 234, 211, 109, 202, 133, 179, 293, 182, 168, 313, 101, 43, 262, 100, 231, 159, 280]
    lifters = [184, 149, 64, 190, 59, 51]
    
    # Sort boxes in descending order
    boxes = [(i, w) for i, w in enumerate(boxes)]
    boxes.sort(key=lambda x: x[1], reverse=True)
    remaining_boxes = boxes.copy()
    solution = []
    
    def can_lift(lifter_indices, weight):
        return sum(lifters[i] for i in lifter_indices) >= weight
    
    while remaining_boxes:
        step = []
        available_lifters = list(range(len(lifters)))
        
        # Try to assign boxes to lifters
        i = 0
        while i < len(remaining_boxes) and available_lifters:
            _, weight = remaining_boxes[i]
            assigned = False
            
            # Try different numbers of lifters
            for r in range(1, min(4, len(available_lifters) + 1)):  # Limit combinations to reduce complexity
                for lifter_combo in itertools.combinations(available_lifters, r):
                    if can_lift(lifter_combo, weight):
                        step.append((weight, list(lifter_combo)))
                        for lifter in lifter_combo:
                            available_lifters.remove(lifter)
                        remaining_boxes.pop(i)
                        assigned = True
                        break
                if assigned:
                    break
            
            if not assigned:
                i += 1
        
        if step:
            solution.append(step)
        else:
            return "Failed to find valid solution"
            
        if len(solution) > 8:
            return "No solution found within 8 steps"
    
    # Format solution
    result = ""
    for i, step in enumerate(solution, 1):
        result += f"Step {i}: {step}\n"
    return result.strip()

print(find_lifting_solution())