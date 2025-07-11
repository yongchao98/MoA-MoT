import collections

# A dictionary to hold the phenotype descriptions for each choice.
phenotype_choices = {
    'A': "Solid black with a metallic blue sheen",
    'B': "Orange with dark blue splotches on the legs",
    'C': "Bright blue with black spots",
    'D': "Green-yellow body with brown or blue legs",
    'E': "Bright red body with blue legs and feet",
    'F': "Solid red",
    'G': "Purple-red body with yellow eyes",
    'H': "Bright blue with white spots across the back",
    'I': "Yellow with black stripes across the back"
}

# The "Isla Colón" morph of Oophaga pumilio is one of the most famous,
# commonly called the "blue jeans" poison dart frog due to its characteristic
# bright red or orange-red body and contrasting blue legs.
correct_choice_letter = 'E'

# Print the final correct choice and its description.
print("The best description for the 'Isla Colón' morph is:")
print(f"Choice {correct_choice_letter}: {phenotype_choices[correct_choice_letter]}")