import numpy as np

def calculate_ideal_ratio():
    """
    Calculates and prints the ideal Ni/Ce molar ratio for catalytic applications.
    
    The ideal ratio for Ni-Ceria nanoparticles in Water Gas Shift and Water Splitting
    is often found where Ni is highly dispersed on the Ceria support. A molar ratio
    of approximately 1:9 (Ni:Ce) is frequently cited as optimal for maximizing
    the active catalytic sites.
    """
    
    # Define the number of atoms based on the commonly cited optimal molar ratio
    ni_atoms = 1
    ce_atoms = 9
    
    # Calculate the decimal ratio
    ratio = ni_atoms / ce_atoms
    
    # Output the explanation and the final equation with each number
    print("The ideal Ni/Ce molar ratio is often found to be approximately 1:9.")
    print(f"This is based on having {ni_atoms} part(s) of Nickel for every {ce_atoms} parts of Cerium.")
    print(f"The calculated decimal ratio is: {ni_atoms} / {ce_atoms} = {ratio:.4f}")

if __name__ == "__main__":
    calculate_ideal_ratio()
