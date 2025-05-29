from itertools import combinations

boxes = [161, 400, 306, 233, 179, 305, 284, 214, 288, 132, 387, 198, 274, 40, 392, 183, 275, 259, 67, 177]
lifters = [122, 116, 51, 165, 95, 150, 179]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Dynamic programming table
dp = {}

def solve(box_index, used_lifters_mask, steps):
    if box_index == len(boxes):
        return steps if steps <= 7 else float('inf')

    if (box_index, used_lifters_mask) in dp:
        return dp[(box_index, used_lifters_mask)]

    min_steps = float('inf')
    for r in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), r):
            if all((used_lifters_mask & (1 << i)) == 0 for i in combo):
                if sum(lifters[i] for i in combo) >= boxes[box_index]:
                    new_mask = used_lifters_mask
                    for i in combo:
                        new_mask |= (1 << i)
                    min_steps = min(min_steps, solve(box_index + 1, new_mask, steps + 1))

    dp[(box_index, used_lifters_mask)] = min_steps
    return min_steps

min_steps = solve(0, 0, 0)

if min_steps > 7:
    print("Unable to lift all boxes within 7 steps.")
else:
    # Reconstruct the solution
    steps = [[] for _ in range(min_steps)]
    box_index = 0
    used_lifters_mask = 0
    current_step = 0

    while box_index < len(boxes):
        for r in range(1, len(lifters) + 1):
            for combo in combinations(range(len(lifters)), r):
                if all((used_lifters_mask & (1 << i)) == 0 for i in combo):
                    if sum(lifters[i] for i in combo) >= boxes[box_index]:
                        new_mask = used_lifters_mask
                        for i in combo:
                            new_mask |= (1 << i)
                        if solve(box_index + 1, new_mask, current_step + 1) == min_steps:
                            steps[current_step].append((boxes[box_index], list(combo)))
                            used_lifters_mask = new_mask
                            box_index += 1
                            break
            else:
                continue
            break
        current_step += 1

    output = "<<<"
    for i, step in enumerate(steps):
        if step:
            output += f"Step {i + 1}: {step}\n"
    output += ">>>"
    print(output)