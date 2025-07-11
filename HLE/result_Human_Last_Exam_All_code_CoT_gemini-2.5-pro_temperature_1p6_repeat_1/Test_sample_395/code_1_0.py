import math

# Define the parameters of the problem
n = 2024  # number of sets
k = 45    # size of each set
l = 1     # size of intersection of any two distinct sets

# Step 1: Calculate the lower bound for the size of the union
# v >= n * k^2 / (n - 1 + k)
v_bound = (n * k**2) / (n - 1 + k)
v_min_theoretical = math.ceil(v_bound)

print(f"n (number of sets) = {n}")
print(f"k (size of each set) = {k}")
print(f"The lower bound from the inequality v >= n*k^2 / (n-1+k) is {v_bound:.3f}")
print(f"As v must be an integer, the minimum possible value is at least {v_min_theoretical}")
print("-" * 20)

# Step 2: Check if v = v_min_theoretical is arithmetically possible
v = v_min_theoretical
# Average replication number r_bar = n*k / v
r_bar = (n * k) / v
print(f"Assuming v = {v}, the average replication number r_bar is {r_bar:.3f}")
print("This suggests the replication numbers r(x) might be integers close to this average.")
print("-" * 20)

# Step 3: Test a possible distribution of replication numbers
# Let's test if r(x) can only take values m and m+1, or m-1 and m+1.
# Based on the calculation in thought, we test r(x) in {44, 46}
r1 = 44
r2 = 46
print(f"We test if a distribution with replication numbers {r1} and {r2} is possible.")

# We have two equations from the double counting relations:
# 1) v_r1 + v_r2 = v
# 2) r1*v_r1 + r2*v_r2 = n*k
S1 = n * k
v_r2 = (S1 - r1 * v) / (r2 - r1)
v_r1 = v - v_r2

if v_r1.is_integer() and v_r1 >= 0 and v_r2.is_integer() and v_r2 >= 0:
    v_r1 = int(v_r1)
    v_r2 = int(v_r2)
    print(f"Solving the system of equations gives:")
    print(f"Number of elements with r(x)={r1} (v_{r1}) = {v_r1}")
    print(f"Number of elements with r(x)={r2} (v_{r2}) = {v_r2}")

    # Step 4: Verify the second moment condition
    # sum(r(x)*(r(x)-1)) = n*(n-1)
    sum_r_r_minus_1 = v_r1 * r1 * (r1 - 1) + v_r2 * r2 * (r2 - 1)
    n_n_minus_1 = n * (n - 1)
    
    print(f"\nChecking the second moment condition:")
    print(f"Calculated sum(r(x)(r(x)-1)) = {v_r1} * {r1} * {r1-1} + {v_r2} * {r2} * {r2-1} = {sum_r_r_minus_1}")
    print(f"Required value n*(n-1) = {n} * {n-1} = {n_n_minus_1}")
    
    if sum_r_r_minus_1 == n_n_minus_1:
        print("\nThe calculated sum matches the required value.")
        print("This distribution is arithmetically consistent.")
        print(f"\nThe smallest possible value of the union is {v}")
    else:
        print("\nThe calculated sum does not match the required value.")
        print("This distribution is not valid. The true minimum might be different or have a different structure for r(x).")
else:
    print(f"No non-negative integer solution found for v_{r1} and v_{r2}.")
    print("This means elements cannot have only replication numbers {r1} and {r2}.")
