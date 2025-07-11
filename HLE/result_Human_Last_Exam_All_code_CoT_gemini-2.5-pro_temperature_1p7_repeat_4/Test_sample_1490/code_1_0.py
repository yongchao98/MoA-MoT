# Step 1: Establish the total number of original curves.
# The prose states the base pattern contains 360 curves. This value allows for integer
# results in the subsequent calculations, resolving contradictions in the poem.
total_original_curves = 360

# Step 2: Calculate the number of curves that need to be redrawn.
# The poem states "one-fourth find new pathways to trace".
curves_to_redraw = total_original_curves * (1/4)

# Step 3: Calculate the breakdown of the newly drawn curves.
# "One-fifth through parabolas sway"
parabolic_curves = int(curves_to_redraw * (1/5))

# "Two-ninths in ellipse arcs reside"
elliptical_curves = int(curves_to_redraw * (2/9))

# "The rest in circles now abide"
circular_curves = int(curves_to_redraw - parabolic_curves - elliptical_curves)

# Step 4: The total number of curves the master must draw is the sum of all new curves.
total_curves_to_draw = parabolic_curves + elliptical_curves + circular_curves

# Final Output: Print the equation showing the sum of the new curves.
print("The master must draw a total of new curves calculated as follows:")
print(f"{parabolic_curves} (parabolic) + {elliptical_curves} (elliptical) + {circular_curves} (circular) = {int(total_curves_to_draw)}")
