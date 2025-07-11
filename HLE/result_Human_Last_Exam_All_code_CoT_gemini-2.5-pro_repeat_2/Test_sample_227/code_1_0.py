import sys

# This script calculates the maximal distance a snail can travel under specific observation conditions.
# The problem is modeled as a system of linear equations derived from a specific arrangement of observers.

# Step 1: Define the number of sub-intervals.
# The 7-minute interval is divided into 21 sub-intervals of 1/3 minute each.
num_intervals = 21

# Step 2: Define an array to hold the distances x_i for each sub-interval.
# We use size num_intervals + 1 for 1-based indexing.
x = [0] * (num_intervals + 1)

# Step 3: Assign values to x_i based on a solution that maximizes the total sum.
# This specific solution satisfies the system of equations derived from the observer arrangement.
# The system of equations is: x_{2k+1} + x_{2k+2} + x_{2k+3} = 1 for k=0..9
x[1] = 1
x[2] = 0
for i in range(3, num_intervals + 1):
    if i % 2 == 1:  # i is odd
        x[i] = 0
    else:  # i is even
        x[i] = 1

# Step 4: Verify that the solution satisfies all 10 equations.
all_ok = True
for k in range(10):
    s = x[2*k+1] + x[2*k+2] + x[2*k+3]
    if s != 1:
        all_ok = False
        print(f"Verification failed for observer k={k}!", file=sys.stderr)
        break

if not all_ok:
     sys.exit("The chosen solution is incorrect.")

# Step 5: Calculate and print the total distance.
# The total distance is the sum of the distances traveled in each sub-interval.
total_distance = sum(x)

# Print the final calculation as an equation with all its terms.
# This shows how the final answer is composed.
print("The maximal possible distance is the sum of the distances from each sub-interval.")
print("Maximal Distance =", " + ".join(map(str, x[1:])), f"= {total_distance}")
