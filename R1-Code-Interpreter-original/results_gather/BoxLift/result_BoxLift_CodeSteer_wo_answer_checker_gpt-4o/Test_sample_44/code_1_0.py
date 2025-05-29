boxes = [116, 246, 296, 369, 78, 275, 77, 383, 71, 155, 352, 94, 63, 168, 350, 79, 59, 252, 88, 278, 188, 383, 240, 308]
lifters = [134, 132, 150, 137, 62, 132]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def assign_lifters(boxes, lifters, step_limit):
    steps = []
    
    def solve(remaining_boxes, current_step):
        if not remaining_boxes:
            return True
        if current_step >= step_limit:
            return False
        
        step = []
        used = [False] * len(lifters)
        
        for box in remaining_boxes[:]:
            lifter_indices = []
            total_capacity = 0
            for i, lifter in enumerate(lifters):
                if not used[i] and total_capacity + lifter <= box:
                    total_capacity += lifter
                    lifter_indices.append(i)
                if total_capacity >= box:
                    break
            
            if total_capacity >= box:
                step.append((box, lifter_indices))
                for idx in lifter_indices:
                    used[idx] = True
                remaining_boxes.remove(box)
        
        if step:
            steps.append(step)
            if solve(remaining_boxes, current_step + 1):
                return True
            steps.pop()
        
        return False
    
    solve(boxes, 0)
    return steps

steps = assign_lifters(boxes, lifters, 9)

# Format the output
output = ""
for i, step in enumerate(steps):
    output += f"Step {i+1}: {step}\n"

print(f"<<<{output.strip()}>>>")