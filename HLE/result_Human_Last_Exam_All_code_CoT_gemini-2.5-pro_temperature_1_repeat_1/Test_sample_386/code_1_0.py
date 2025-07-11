import re

def solve_cell_culture_problem():
    """
    This function analyzes the correct statement about isolating a fibroblast cell line
    and extracts the relevant numerical data for the culture medium composition.
    """
    # The correct statement is C.
    correct_statement = "Debrided epithelium and endothelium induced proliferation of the stromal cells to myofibroblasts in the medium containing 10% of the FBS and 1% antibiotic, and they adhered to the bottom of the flask."

    # Use regular expressions to find all numbers in the string.
    numbers = re.findall(r'\d+', correct_statement)

    # Convert the found strings to integers.
    fbs_percentage = int(numbers[0])
    antibiotic_percentage = int(numbers[1])

    # Print the final equation representing the medium composition.
    print("Based on the correct procedure described in option C, the culture medium composition is:")
    print(f"Final Medium = {fbs_percentage}% FBS + {antibiotic_percentage}% Antibiotic")

solve_cell_culture_problem()