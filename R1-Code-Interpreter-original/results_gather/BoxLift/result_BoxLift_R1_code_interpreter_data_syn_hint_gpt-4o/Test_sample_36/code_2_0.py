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
    
    def construct_step(remaining_boxes, available_lifters):
        step = []
        for box in remaining_boxes[:]:
            for r in range(1, len(available_lifters) + 1):
                for lifter_combination in combinations(available_lifters, r):
                    if can_lift(box, lifter_combination):
                        step.append((box, list(lifter_combination)))
                        for i in lifter_combination:
                            available_lifters.remove(i)
                        remaining_boxes.remove(box)
                        break
                else:
                    continue
                break
        return step
    
    while boxes:
        available_lifters = lifter_indices[:]
        step = construct_step(boxes, available_lifters)
        if not step:
            return "Not possible to lift all boxes in 6 steps"
        steps.append(step)
        if len(steps) > 6:
            return "Not possible to lift all boxes in 6 steps"
    
    return steps

# Get the steps
steps = assign_lifters_to_boxes(boxes, lifters)

# Print the steps
print(steps)