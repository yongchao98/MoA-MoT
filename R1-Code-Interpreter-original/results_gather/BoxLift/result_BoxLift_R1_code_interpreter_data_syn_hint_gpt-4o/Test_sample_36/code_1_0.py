from itertools import combinations

# Box weights and lifter capacities
boxes = [85, 162, 147, 83, 142, 96, 200, 172, 151, 77, 59, 39]
lifters = [53, 54, 76, 48, 97]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Function to find the minimum number of steps
def assign_lifters_to_boxes(boxes, lifters):
    steps = []
    lifter_indices = list(range(len(lifters)))
    
    def can_lift(box, lifter_combination):
        return sum(lifters[i] for i in lifter_combination) >= box
    
    def backtrack(remaining_boxes, current_step):
        if not remaining_boxes:
            steps.append(current_step)
            return True
        if len(steps) >= 6:
            return False
        
        for box in remaining_boxes:
            for r in range(1, len(lifters) + 1):
                for lifter_combination in combinations(lifter_indices, r):
                    if can_lift(box, lifter_combination):
                        new_step = current_step + [(box, list(lifter_combination))]
                        new_remaining_boxes = remaining_boxes[:]
                        new_remaining_boxes.remove(box)
                        if backtrack(new_remaining_boxes, new_step):
                            return True
        return False
    
    if backtrack(boxes, []):
        return steps
    else:
        return "Not possible to lift all boxes in 6 steps"

# Get the steps
steps = assign_lifters_to_boxes(boxes, lifters)

# Print the steps
print(steps)