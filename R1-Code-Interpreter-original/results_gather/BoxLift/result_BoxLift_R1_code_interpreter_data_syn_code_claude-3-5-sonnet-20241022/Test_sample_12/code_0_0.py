from itertools import combinations

def can_lift(lifters, weight):
    return sum(lifters) >= weight

def find_minimal_lifters(available_lifters, weight):
    for r in range(1, len(available_lifters) + 1):
        for combo in combinations(available_lifters, r):
            if sum(combo) >= weight:
                return list(combo)
    return None

def solve_box_lifting():
    boxes = [50, 219, 120, 245, 204, 172, 274, 225, 94, 270, 158, 253, 174, 40, 252, 124]
    lifter_capacities = [54, 61, 106, 60, 63]
    
    # Sort boxes in descending order
    boxes = sorted(boxes, reverse=True)
    
    steps = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes:
        step = []
        used_lifters = set()
        
        # Try to assign as many boxes as possible in current step
        for box in remaining_boxes[:]:
            available_lifters = [(i, cap) for i, cap in enumerate(lifter_capacities) 
                               if i not in used_lifters]
            
            if not available_lifters:
                break
                
            # Get lifter indices and their capacities
            lifter_indices = [x[0] for x in available_lifters]
            lifter_caps = [x[1] for x in available_lifters]
            
            # Find minimal combination of lifters for this box
            lifters_needed = find_minimal_lifters(lifter_caps, box)
            
            if lifters_needed:
                # Map capacities back to indices
                selected_lifters = []
                for cap in lifters_needed:
                    idx = lifter_caps.index(cap)
                    selected_lifters.append(lifter_indices[idx])
                    used_lifters.add(lifter_indices[idx])
                    lifter_caps.pop(idx)
                    lifter_indices.pop(idx)
                
                step.append((box, sorted(selected_lifters)))
                remaining_boxes.remove(box)
        
        steps.append(step)
    
    # Format output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    print(result.strip())

solve_box_lifting()