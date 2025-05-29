boxes = [76, 97, 14, 30, 36, 83, 45, 46]
lifters = [69, 59, 68, 50]

# Sort boxes in descending order
boxes.sort(reverse=True)

def find_combination(box, lifters, used, current_combination, start):
    if sum(current_combination) >= box:
        return current_combination
    for i in range(start, len(lifters)):
        if not used[i]:
            used[i] = True
            result = find_combination(box, lifters, used, current_combination + [lifters[i]], i + 1)
            if result:
                return result
            used[i] = False
    return None

steps = []

for step in range(3):
    step_lift = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:
        lifter_indices = []
        combination = find_combination(box, lifters, used_lifters, [], 0)
        if combination:
            for weight in combination:
                index = lifters.index(weight)
                lifter_indices.append(index)
                used_lifters[index] = True
            step_lift.append((box, lifter_indices))
            boxes.remove(box)
    
    if step_lift:
        steps.append(f"Step {step + 1}: {step_lift}")
    if not boxes:
        break

result = "\n".join(steps)
print(f"<<<{result}>>>")