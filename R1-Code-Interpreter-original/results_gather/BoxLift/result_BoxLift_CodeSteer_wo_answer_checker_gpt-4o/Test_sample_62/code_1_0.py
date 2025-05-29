from pulp import LpProblem, LpMinimize, LpVariable, lpSum, LpBinary, value

# Data
boxes = [91, 207, 152, 47, 209, 54, 251, 176, 194, 221, 152, 141, 128, 159, 57, 184]
lifters = [149, 131, 113, 109, 124]
max_steps = 6

# Problem
prob = LpProblem("BoxLift", LpMinimize)

# Variables
x = LpVariable.dicts("x", ((i, j, k) for i in range(len(boxes)) for j in range(len(lifters)) for k in range(max_steps)), 0, 1, LpBinary)

# Objective: Minimize the number of steps
prob += lpSum(x[i, j, k] for i in range(len(boxes)) for j in range(len(lifters)) for k in range(max_steps))

# Constraints
# Each box must be lifted exactly once
for i in range(len(boxes)):
    prob += lpSum(x[i, j, k] for j in range(len(lifters)) for k in range(max_steps)) == 1

# Each lifter can be used only once per step
for j in range(len(lifters)):
    for k in range(max_steps):
        prob += lpSum(x[i, j, k] for i in range(len(boxes))) <= 1

# Total weight lifted in a step must not exceed lifters' capacity
for k in range(max_steps):
    for i in range(len(boxes)):
        prob += lpSum(x[i, j, k] * lifters[j] for j in range(len(lifters))) >= boxes[i] * lpSum(x[i, j, k] for j in range(len(lifters)))

# Solve
prob.solve()

# Output
steps = [[] for _ in range(max_steps)]
for i in range(len(boxes)):
    for j in range(len(lifters)):
        for k in range(max_steps):
            if value(x[i, j, k]) == 1:
                steps[k].append((boxes[i], [j]))

output = []
for i, step in enumerate(steps):
    if step:
        output.append(f"Step {i+1}: {step}")

print("<<<" + "\n".join(output) + ">>>")