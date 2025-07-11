# Problem parameters
total_duration = 7  # minutes

# Optimal gadget parameters
gadget_duration = 1.5  # minutes
gadget_distance = 2  # meters

# Step 1: Determine how many times the gadget can be used to tile the total duration.
num_gadgets = int(total_duration / gadget_duration)

# Step 2: Calculate the duration covered by the gadgets and the distance traveled.
duration_covered_by_gadgets = num_gadgets * gadget_duration
distance_from_gadgets = num_gadgets * gadget_distance

# Step 3: Calculate the remaining duration.
remaining_duration = total_duration - duration_covered_by_gadgets

# Step 4: The distance traveled in the remaining time is equal to the time elapsed,
# as it's a single block covered by a single observer.
distance_from_remainder = remaining_duration

# Step 5: Calculate the total maximal distance.
total_distance = distance_from_gadgets + distance_from_remainder

# Step 6: Print the explanation and the final equation.
print("The maximal distance can be found by tiling the 7-minute interval with an optimal 'gadget'.")
print("An optimal gadget covers 1.5 minutes and allows the snail to travel 2 meters.")
print(f"The 7-minute interval can be broken down into {num_gadgets} such gadgets and a remainder.")
print(f"Duration covered by gadgets: {num_gadgets} * {gadget_duration} = {duration_covered_by_gadgets} minutes.")
print(f"Distance from these gadgets: {num_gadgets} * {gadget_distance} = {distance_from_gadgets} meters.")
print(f"Remaining duration: {total_duration} - {duration_covered_by_gadgets} = {remaining_duration} minute.")
print(f"Distance covered in the remaining minute: {distance_from_remainder} meter.")
print("\nThe final calculation for the maximal distance is:")

# Output the equation with each number explicitly printed
print(f"{num_gadgets} * {gadget_distance} + {int(distance_from_remainder)} = {int(total_distance)}")
