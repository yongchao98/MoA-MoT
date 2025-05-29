from itertools import combinations

boxes = [44, 178, 94, 176, 159, 199, 197, 133, 173, 186, 196, 80]
lifters = [61, 45, 66, 38, 52]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Function to find the minimum number of steps
def find_min_steps(boxes, lifters):
    steps = []
    while boxes:
        step = []
        used_lifters = [False] * len(lifters)
        
        for box in boxes[:]:
            for r in range(1, len(lifters) + 1):
                for combo in combinations(enumerate(lifters), r):
                    indices, capacities = zip(*combo)
                    if not any(used_lifters[i] for i in indices) and sum(capacities) >= box:
                        step.append((box, list(indices)))
                        for i in indices:
                            used_lifters[i] = True
                        boxes.remove(box)
                        break
                else:
                    continue
                break
        
        steps.append(step)
    
    return steps

# Get the steps
steps = find_min_steps(boxes, lifters)

# Print the steps
output = []
for i, step in enumerate(steps):
    output.append(f"Step {i+1}: {step}")

print("<<<" + "\n".join(output) + ">>>")