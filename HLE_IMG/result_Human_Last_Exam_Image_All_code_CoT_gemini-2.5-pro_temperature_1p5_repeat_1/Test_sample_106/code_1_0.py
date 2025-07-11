# Based on the analysis of the Wagner-Meerwein rearrangement mechanism,
# we can determine the identity of the substituents at each numbered position.
substituents = {
    1: 'CH3',
    2: 'H',
    3: 'CH3',
    4: 'CH3',
    5: 'H'
}

# Create a list of strings in the format "position = substituent"
# e.g., ["1 = CH3", "2 = H", ...]
answer_parts = [f"{position} = {substituent}" for position, substituent in substituents.items()]

# Join the parts with a comma and a space for the final output string
final_answer_string = ", ".join(answer_parts)

print(final_answer_string)