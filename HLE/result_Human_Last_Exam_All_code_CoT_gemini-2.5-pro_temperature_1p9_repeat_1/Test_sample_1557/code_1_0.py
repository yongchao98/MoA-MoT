# The task is to identify the correct statements about Duke of Wellington's career.
# Based on a historical analysis, the correct statements have been identified.
# This script will print the numbers of the correct statements, sorted and separated by a comma.

# List of the numbers corresponding to the correct statements.
correct_options = [1, 6, 8]

# Sort the list to ensure the output is in the correct order.
correct_options.sort()

# Convert the numbers to strings and join them with a comma and space for the final output.
# The question requests the output of "each number in the final equation", 
# which is interpreted as printing the final list of numbers clearly.
final_answer = ", ".join(map(str, correct_options))

print("The numbers of the correct statements are:")
print(final_answer)

print("\nFor direct submission:")
# The format requested is <<<answer content>>>
print(f"<<<{final_answer}>>>")