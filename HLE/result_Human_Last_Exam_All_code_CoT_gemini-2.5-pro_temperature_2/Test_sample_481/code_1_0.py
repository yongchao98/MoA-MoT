# This script identifies and displays the different numbers imprinted on
# Teva brand fluoxetine capsules for 10mg and 20mg dosages,
# based on publicly available drug identification data.

# Assign the identified numbers to variables.
imprint_num_10mg = 3105
imprint_num_20mg = 3109

print("The following are the different numbers imprinted on the capsules:")
print(f"Teva brand fluoxetine 10mg capsule has the number: {imprint_num_10mg}")
print(f"Teva brand fluoxetine 20mg capsule has the number: {imprint_num_20mg}")

# Final Answer formatting
# The question asks for the "different numbers".
final_answer_string = f"{imprint_num_10mg} and {imprint_num_20mg}"
print(f"\nIn summary, the two different numbers are {final_answer_string}.")

# The final answer for the wrapper.
print(f"<<<{final_answer_string}>>>")