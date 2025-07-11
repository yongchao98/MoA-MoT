# Step 1: Define the initial values from the problem description.
initial_curves = 360
unchanged_curves = 90

print(f"The Rangoli pattern initially has {initial_curves} curves.")
print(f"After the disruption, {unchanged_curves} curves remain unchanged.")
print("-" * 20)

# Step 2: Calculate the number of curves that were disturbed and need to be replaced.
# The poem states "Three-eighths of curves lost shape serene".
disturbed_curves = (3/8) * initial_curves

# We must have a whole number of curves.
disturbed_curves = int(disturbed_curves)
print(f"Number of curves that need to be replaced: (3/8) * {initial_curves} = {disturbed_curves}")
print("-" * 20)

# Step 3: Calculate the number of each new type of curve.
# The poem provides the breakdown for these disturbed curves.
# "One-fifth through parabolas sway"
parabolic_curves = int((1/5) * disturbed_curves)
print(f"Number of new parabolic curves: (1/5) * {disturbed_curves} = {parabolic_curves}")

# "Two-ninths in ellipse arcs reside"
elliptical_curves = int((2/9) * disturbed_curves)
print(f"Number of new elliptical curves: (2/9) * {disturbed_curves} = {elliptical_curves}")

# "The rest in circles now abide"
circular_curves = disturbed_curves - parabolic_curves - elliptical_curves
print(f"Number of new circular curves: {disturbed_curves} - {parabolic_curves} - {elliptical_curves} = {circular_curves}")
print("-" * 20)

# Step 4: Calculate the total number of curves in the final, restored pattern.
# This is the sum of the unchanged curves and all the new curves drawn by the master.
total_final_curves = unchanged_curves + parabolic_curves + elliptical_curves + circular_curves

# Step 5: Print the final equation and the answer.
print("To restore the pattern, the total number of curves is the sum of the unchanged and all the new curves.")
print("Final Equation:")
print(f"{unchanged_curves} (unchanged) + {parabolic_curves} (parabolic) + {elliptical_curves} (elliptical) + {circular_curves} (circular) = {total_final_curves}")

print("\nFinal Answer:")
print(f"The master must draw a total of {total_final_curves} curves to restore the pattern.")

<<<225>>>