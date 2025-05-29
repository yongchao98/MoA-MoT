from pulp import LpProblem, LpMinimize, LpVariable, lpSum, LpBinary, value

boxes = [289, 375, 107, 145, 257, 48, 141, 83, 136, 368, 59, 133, 186, 266, 353, 73, 66, 210, 247, 79, 342, 318, 337, 162]
lifters = [75, 189, 162, 181, 118, 194, 137]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Number of boxes and lifters
num_boxes = len(boxes)
num_lifters = len(lifters)
max_steps = 7

# Create the problem
prob = LpProblem("BoxLift", LpMinimize)

# Create variables
x = LpVariable.dicts("x", (range(num_boxes), range(num_lifters), range(max_steps)), 0, 1, LpBinary)

# Objective function: minimize the number of steps
prob += lpSum(x[i][j][k] for i in range(num_boxes) for j in range(num_lifters) for k in range(max_steps))

# Constraints
# Each box must be lifted exactly once
for i in range(num_boxes):
    prob += lpSum(x[i][j][k] for j in range(num_lifters) for k in range(max_steps)) == 1

# The total weight lifted by the selected lifters for a box must be at least the weight of the box
for i in range(num_boxes):
    for k in range(max_steps):
        prob += lpSum(lifters[j] * x[i][j][k] for j in range(num_lifters)) >= boxes[i]

# Each lifter can be used at most once per step
for j in range(num_lifters):
    for k in range(max_steps):
        prob += lpSum(x[i][j][k] for i in range(num_boxes)) <= 1

# Solve the problem
prob.solve()

# Extract the solution
steps = [[] for _ in range(max_steps)]
for i in range(num_boxes):
    for j in range(num_lifters):
        for k in range(max_steps):
            if value(x[i][j][k]) == 1:
                steps[k].append((boxes[i], [j]))

# Filter out empty steps and print the result
output = []
for i, step in enumerate(steps):
    if step:
        output.append(f"Step {i+1}: {step}")

print("<<<" + "\n".join(output) + ">>>")