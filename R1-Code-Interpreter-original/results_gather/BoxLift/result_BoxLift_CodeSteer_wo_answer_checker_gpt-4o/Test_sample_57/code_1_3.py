from itertools import combinations

boxes = [34, 92, 289, 35, 140, 174, 244, 136, 80, 252, 136, 181, 92, 235, 61, 286]
lifters = [61, 60, 150, 96, 58, 91]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Memoization dictionary
memo = {}

def can_lift(box, lifters):
    for r in range(1, len(lifters) + 1):
        for combo in combinations(lifters, r):
            if sum(combo) >= box:
                return list(combo)
    return None

def backtrack(remaining_boxes, available_lifters, steps):
    if not remaining_boxes:
        return steps
    
    if len(steps) >= 6:
        return None
    
    state = (tuple(remaining_boxes), tuple(available_lifters))
    if state in memo:
        return memo[state]
    
    for i, box in enumerate(remaining_boxes):
        lifter_combo = can_lift(box, available_lifters)
        if lifter_combo:
            new_lifters = available_lifters[:]
            for lifter in lifter_combo:
                new_lifters.remove(lifter)
            
            new_steps = steps[:]
            new_steps.append((box, [lifters.index(l) for l in lifter_combo]))
            
            result = backtrack(remaining_boxes[:i] + remaining_boxes[i+1:], lifters, new_steps)
            if result:
                memo[state] = result
                return result
    
    memo[state] = None
    return None

# Start backtracking
solution = backtrack(boxes, lifters, [])

# Print the solution
if solution:
    output = []
    for i, step in enumerate(solution):
        output.append(f"Step {i+1}: {step}")
    print("<<<" + "\n".join(output) + ">>>")
else:
    print("No solution found within 6 steps.")