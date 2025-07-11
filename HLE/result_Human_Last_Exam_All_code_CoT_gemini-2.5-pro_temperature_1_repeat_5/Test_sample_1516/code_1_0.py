import math

# Define the given parameters
# r0 is the initial radius of the first row in meters.
r0 = 3
# h_seated is the height of a seated person in meters.
h_seated = 1.0
# h_standing is the height of a standing person in meters.
h_standing = 1.5

print("Step 1: Identify the critical visibility constraint.")
print("The explicit constraint of seeing a speaker in the front row is always satisfied for any positive K.")
print("Therefore, we use a more general and restrictive constraint: every person must be able to see over the head of the person in front of them towards a central point (0,0).\n")

print("Step 2: Formulate the constraint mathematically.")
print("The head of a person in a row at radius 'r' is at height h(r) = (r^2 / K) + h_seated.")
print("The slope of the sightline to the center (0,0) is [h(r) + h_seated] / r.")
print("For visibility, this slope must increase with r. So, the function f(r) = (r^2/K + 1)/r must be increasing.")
print("f(r) = r/K + 1/r\n")


print("Step 3: Solve the inequality for K.")
print("For f(r) to be increasing, its derivative f'(r) must be greater than 0.")
print("f'(r) = 1/K - 1/r^2 > 0")
print("This simplifies to 1/K > 1/r^2, which means K < r^2.\n")

print("Step 4: Apply the most restrictive case.")
print("This condition K < r^2 must hold for all radii 'r' where people are seated.")
print("The most restrictive case is for the smallest radius, which is the first row.")
print(f"The radius of the first row (r0) is {r0} m.\n")

print("Step 5: Calculate the final value for K.")
r0_squared = r0**2
print(f"The final inequality is K < r0^2, which is:")
print(f"K < {r0}^2")
print(f"K < {r0_squared}")

# Since K must be an integer, the maximum value is the integer just below r0_squared.
max_K = math.floor(r0_squared - 1)
print(f"As K must be an integer, the maximum value K can take is {max_K}.")

print("\nFinal Answer:")
print(f"The final equation is K < {r0}^2, which means K < {r0_squared}. The maximum integer value for K is {max_K}.")
print(f"<<<{max_K}>>>")
