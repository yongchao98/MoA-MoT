# This script calculates the time for an egg to fall from a skyscraper.
# The calculation uses simple integers and fractions for a young learner.

print("Let's figure out how long it takes for the egg to reach the ground!")
print("First, we need to know how tall the skyscraper is.\n")

# --- Part 1: Finding the height of the skyscraper ---
# Knowns for height calculation
distance = 240
# We use the approximation tan(22.5 degrees) ≈ 5/12
tan_approx_num = 5
tan_approx_den = 12

print(f"We are standing {distance} meters away. We can find the height 'h' with this calculation:")
print(f"h = {distance} * {tan_approx_num} / {tan_approx_den}")

# To make it easier, do the division first
h_intermediate = distance / tan_approx_den
print(f"Let's calculate {distance} / {tan_approx_den} first, which is {int(h_intermediate)}.")

height = h_intermediate * tan_approx_num
print(f"Then, we multiply by {tan_approx_num}. The height is {int(height)} meters.\n")


# --- Part 2: Finding the time of the fall ---
print("Now that we know the height is 100 meters, we can find the fall time 't'.")

# Approximation for gravity, g ≈ 10 m/s^2
g_approx = 10
print(f"The formula for 'time squared' (t*t) is 2 * height / g.")
print(f"We'll use a simple number for g: {g_approx}.")

t_squared = 2 * height / g_approx
print(f"So, t*t = 2 * {int(height)} / {g_approx} = {int(t_squared)}\n")


# --- Part 3: Finding the final time ---
print("This means the time 't' is the square root of 20.")
print("We can simplify this! 20 is 4 times 5, so the square root of 20 is 2 times the square root of 5.\n")

# Approximation for sqrt(5) ≈ 9/4
sqrt5_approx_num = 9
sqrt5_approx_den = 4
print(f"A good and simple fraction for the square root of 5 is {sqrt5_approx_num}/{sqrt5_approx_den}.")

print("\nHere is the final calculation for the time, which is our answer:")
# The final equation for time is t ≈ 2 * (9/4)
print(f"Time = 2 * {sqrt5_approx_num} / {sqrt5_approx_den}")

# Final result
final_time_numerator = 2 * sqrt5_approx_num
final_time_denominator = sqrt5_approx_den
final_time = final_time_numerator / final_time_denominator
print(f"So the total time is {final_time} seconds.")