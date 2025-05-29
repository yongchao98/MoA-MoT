import heapq

boxes = [88, 128, 39, 30, 80, 60, 162, 154, 191, 150, 178, 32]
lifters = [95, 111, 115, 35]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Improved heuristic function
def heuristic(remaining_boxes, lifters):
    total_weight = sum(remaining_boxes)
    max_capacity = sum(lifters)
    return (total_weight + max_capacity - 1) // max_capacity

# Branch and bound approach with refined heuristic
def branch_and_bound(boxes, lifters):
    # Priority queue for exploring paths
    pq = []
    heapq.heappush(pq, (0, [], boxes, [False] * len(lifters)))
    
    while pq:
        steps_count, steps, remaining_boxes, used_lifters = heapq.heappop(pq)
        
        if not remaining_boxes:
            return steps
        
        if steps_count >= 5:
            continue
        
        for i in range(len(remaining_boxes)):
            box = remaining_boxes[i]
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
                new_boxes = remaining_boxes[:i] + remaining_boxes[i+1:]
                estimated_steps = steps_count + 1 + heuristic(new_boxes, lifters)
                
                if estimated_steps <= 5:
                    heapq.heappush(pq, (steps_count + 1, new_steps, new_boxes, [False] * len(lifters)))
            
            for j in lifter_indices:
                used_lifters[j] = False
    
    return None

# Find the solution
solution = branch_and_bound(boxes, lifters)

# Print the solution
if solution:
    output = []
    for i, step in enumerate(solution):
        output.append(f"Step {i + 1}: {step}")
    print("<<<" + "\n".join(output) + ">>>")
else:
    print("No solution found within 5 steps.")