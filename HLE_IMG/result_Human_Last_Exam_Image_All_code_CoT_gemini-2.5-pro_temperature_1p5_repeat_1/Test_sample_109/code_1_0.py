import math

# Step 1: Define the knowns and derived dimensions.
r_y = 1  # Radius of a yellow circle
R_w = 4  # Radius of a white circle
W = 34   # Width of the target
H = 24   # Height of the target
N_y = 5  # Number of yellow circles
shots = 10000

# Step 2: Calculate the total area of the target.
total_area = W * H
print(f"The total area of the target is the width multiplied by the height.")
print(f"Area = {W} cm * {H} cm = {total_area} cm^2")
print("-" * 20)

# Step 3: Calculate the total area of the yellow circles.
area_one_yellow = math.pi * (r_y ** 2)
total_yellow_area = N_y * area_one_yellow
print(f"There are {N_y} yellow circles, each with a radius of {r_y} cm.")
print(f"The area of one yellow circle is pi * r^2 = pi * {r_y}^2 = {area_one_yellow:.4f} cm^2.")
print(f"The total area of all yellow circles is {N_y} * pi = {total_yellow_area:.4f} cm^2.")
print("-" * 20)

# Step 4: Calculate the probability of hitting a yellow circle.
probability = total_yellow_area / total_area
print("The probability of hitting a yellow circle is the ratio of the yellow area to the total area.")
print(f"P(hit) = ({N_y} * pi) / {total_area} = {probability:.6f}")
print("-" * 20)

# Step 5: Calculate the expected number of hits in 10000 shots.
expected_hits = shots * probability
print("The expected number of hits is the total number of shots multiplied by the probability.")
print(f"Expected Hits = {shots} * (({N_y} * pi) / {total_area})")
print(f"Expected Hits = {expected_hits:.4f}")

# Final Answer as a rounded integer, as "number of times" suggests.
final_answer = round(expected_hits)
print(f"\nRounding to the nearest whole number, we expect to hit a yellow circle approximately {final_answer} times.")
