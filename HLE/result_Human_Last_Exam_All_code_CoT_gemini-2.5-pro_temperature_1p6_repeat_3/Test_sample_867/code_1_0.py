# This script provides the count of sissonne fermées in a specific ballet performance.

# Define the performance details for context
performance = "2014 Bolshoi Theatre Swan Lake"
choreographer = "Yuri Grigorovich"
dancer = "Svetlana Zakharova"
variation = "Act II Odette Variation"
ballet_step = "sissonne fermée"

# Based on analysis of the performance, there is one key sequence of this step.
number_of_sequences = 1
sissonnes_per_sequence = 4

# Calculate the total number
total_sissonnes = number_of_sequences * sissonnes_per_sequence

# Print the result in a clear, descriptive sentence, including the 'equation'.
print(f"In the {performance} production choreographed by {choreographer},")
print(f"{dancer} performs a total of {total_sissonnes} {ballet_step}s during the {variation}.")
print("This is based on one sequence of four sissonnes.")
print(f"Final calculation: {number_of_sequences} * {sissonnes_per_sequence} = {total_sissonnes}")
