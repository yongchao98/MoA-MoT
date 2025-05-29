from pulp import LpProblem, LpMinimize, LpVariable, lpSum, LpBinary, value, LpStatus

# Data
boxes = [57, 112, 121, 49, 166, 143, 169, 43, 90, 66, 23, 143]
lifters = [103, 49, 95, 97, 98]
max_steps = 4

# Problem
prob = LpProblem("BoxLift", LpMinimize)

# Variables
x = LpVariable.dicts("x", ((i, j, k) for i in range(len(boxes)) for j in range(len(lifters)) for k in range(max_steps)), cat=LpBinary)

# Objective: Minimize the number of steps
prob += lpSum(x[i, j, k] for i in range(len(boxes)) for j in range(len(lifters)) for k in range(max_steps))

# Constraints
# Each box must be lifted exactly once
for i in range(len(boxes)):
    prob += lpSum(x[i, j, k] for j in range(len(lifters)) for k in range(max_steps)) == 1

# Each lifter can be used at most once per step
for j in range(len(lifters)):
    for k in range(max_steps):
        prob += lpSum(x[i, j, k] for i in range(len(boxes))) <= 1

# Total weight lifted by assigned lifters must meet or exceed the weight of each box
for i in range(len(boxes)):
    for k in range(max_steps):
        prob += lpSum(x[i, j, k] * lifters[j] for j in range(len(lifters))) >= boxes[i]

# Solve the problem
prob.solve()

# Check if a solution was found
if LpStatus[prob.status] == 'Optimal':
    # Extract the solution
    steps = [[] for _ in range(max_steps)]
    for i in range(len(boxes)):
        for j in range(len(lifters)):
            for k in range(max_steps):
                if value(x[i, j, k]) == 1:
                    steps[k].append((boxes[i], [j]))

    # Filter out empty steps
    steps = [step for step in steps if step]

    # Print the steps
    output = []
    for i, step in enumerate(steps):
        output.append(f"Step {i+1}: {step}")

    print("<<<" + "\n".join(output) + ">>>")
else:
    print("No solution found within 4 steps.")