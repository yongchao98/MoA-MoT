boxes = [88, 128, 39, 30, 80, 60, 162, 154, 191, 150, 178, 32]
lifters = [95, 111, 115, 35]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Function to find a valid solution using backtracking
def find_solution(boxes, lifters, step, steps, used_lifters):
    if not boxes:
        return steps if len(steps) <= 5 else None
    
    if step >= 5:
        return None
    
    for i in range(len(boxes)):
        box = boxes[i]
        lifter_indices = []
        total_capacity = 0
        
        for j, capacity in enumerate(lifters):
            if not used_lifters[j] and total_capacity < box:
                lifter_indices.append(j)
                total_capacity += capacity
                used_lifters[j] = True
            
            if total_capacity >= box:
                break
        
        if total_capacity >= box:
            new_steps = steps[:]
            new_steps.append((box, lifter_indices))
            new_boxes = boxes[:i] + boxes[i+1:]
            result = find_solution(new_boxes, lifters, step + 1, new_steps, [False] * len(lifters))
            if result:
                return result
        
        for j in lifter_indices:
            used_lifters[j] = False
    
    return None

# Initialize used lifters
used_lifters = [False] * len(lifters)

# Find the solution
solution = find_solution(boxes, lifters, 0, [], used_lifters)

# Print the solution
if solution:
    output = []
    for i, step in enumerate(solution):
        output.append(f"Step {i + 1}: {step}")
    print("<<<" + "\n".join(output) + ">>>")
else:
    print("No solution found within 5 steps.")