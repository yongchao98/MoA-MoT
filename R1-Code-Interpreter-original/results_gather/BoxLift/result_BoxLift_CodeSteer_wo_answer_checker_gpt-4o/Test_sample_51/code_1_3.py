boxes = [142, 109, 62, 158, 37, 50, 180, 128, 61, 123, 44, 51]
lifters = [93, 71, 82, 91]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def find_steps(boxes, lifters, max_steps):
    steps = []
    used_lifters = [False] * len(lifters)

    def backtrack(box_index, current_step):
        if box_index == len(boxes):
            steps.append(current_step[:])
            return True

        box = boxes[box_index]
        for i in range(len(lifters)):
            if not used_lifters[i] and lifters[i] >= box:
                used_lifters[i] = True
                current_step.append((box, [i]))
                if backtrack(box_index + 1, current_step):
                    return True
                current_step.pop()
                used_lifters[i] = False

        for i in range(len(lifters)):
            if not used_lifters[i]:
                used_lifters[i] = True
                current_step.append((box, [i]))
                if backtrack(box_index + 1, current_step):
                    return True
                current_step.pop()
                used_lifters[i] = False

        return False

    for _ in range(max_steps):
        if backtrack(0, []):
            break

    return steps

steps = find_steps(boxes, lifters, 5)

# Print the steps
output = ""
for i, step in enumerate(steps):
    output += f"Step {i + 1}: {step}\n"

print(f"<<<{output.strip()}>>>")