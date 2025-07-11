import re

# The text of the correct answer choice
correct_answer_text = "Debrided epithelium and endothelium induced proliferation of the stromal cells to myofibroblasts in the medium containing 10% of the FBS and 1% antibiotic, and they adhered to the bottom of the flask."

# Find all numbers that are followed by a '%' sign
percentages = re.findall(r'(\d+)%', correct_answer_text)

# Check if we found the expected numbers
if len(percentages) == 2:
    fbs_percentage = percentages[0]
    antibiotic_percentage = percentages[1]
    
    # Print the "final equation" for the medium supplement
    print("Culture Medium Supplement Formula:")
    print(f"Final Concentration of FBS = {fbs_percentage}%")
    print(f"Final Concentration of Antibiotic = {antibiotic_percentage}%")
else:
    print("Could not extract the expected percentages from the text.")
