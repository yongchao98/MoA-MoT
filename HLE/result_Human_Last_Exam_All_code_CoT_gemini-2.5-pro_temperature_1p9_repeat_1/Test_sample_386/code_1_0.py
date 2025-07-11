import re

# The text from the correct answer choice (C)
correct_procedure_text = "Debrided epithelium and endothelium induced proliferation of the stromal cells to myofibroblasts in the medium containing 10% of the FBS and 1% antibiotic, and they adhered to the bottom of the flask."

# Find all numbers in the text
numbers = re.findall(r'\d+', correct_procedure_text)

# Convert found numbers from strings to integers
fbs_percentage = int(numbers[0])
antibiotic_percentage = int(numbers[1])

# Calculate the total percentage of the specified additives
total_additives = fbs_percentage + antibiotic_percentage

# Display the information in the requested equation format
print("Analysis of medium composition from the correct procedure:")
print(f"The equation for the primary medium additives is:")
print(f"{fbs_percentage}% (FBS) + {antibiotic_percentage}% (Antibiotic) = {total_additives}% (Total Additives)")