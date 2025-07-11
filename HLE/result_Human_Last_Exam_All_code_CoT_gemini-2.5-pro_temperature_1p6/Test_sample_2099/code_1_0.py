import numpy as np
from mendeleev import element
import re

def solve_saber_duel_riddle():
    """
    Solves the problem by deciphering the elements based on the problem's context,
    calculating the specified graph energy for each resulting molecule, identifying
    the molecule with the lowest energy, and computing the final product.
    """
    # Step 1: Decipher the elements X1-X14 by mapping them to the 14 elements of Period 4,
    # a choice motivated by the "weapon preparation" theme.
    period_4_symbols = ['K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge']
    x_map = {f"X{i+1}": symbol for i, symbol in enumerate(period_4_symbols)}
    
    # Pre-fetch atomic mass data for all required elements.
    atomic_masses = {el.symbol: el.atomic_weight for el in element(period_4_symbols)}

    # Step 2: Define the molecular formulas for Y1-Y10 from the problem statement.
    y_formulas_str = {
        "Y1": "X1X2X3X4X5X6",
        "Y2": "X8X9X10X11X5X6X12X8X9X8X12X13",
        "Y3": "X11X3X4X14X3X4X5X6",
        "Y4": "X12X4X13X5X6X3",
        "Y5": "X8X9X10X11X9X14X5X6X3",
        "Y6": "X1X10X5X1X9X4X3",
        "Y7": "X8X9X10X11X12X4X5X6",
        "Y8": "X10X2X5X13X9X4X12X4X3",
        "Y9": "X9X14X5X11X3X4X14X3X4X3",
        "Y10": "X1X12X1X3X10X12X13X12X4X3",
    }

    # Convert the formula strings into lists of element symbols for each molecule.
    y_molecules = {}
    for name, formula_str in y_formulas_str.items():
        # Use regex to correctly parse variables like X1, X10, etc.
        x_vars = re.findall(r'X\d+', formula_str)
        y_molecules[name] = [x_map[var] for var in x_vars]

    # Step 3 & 4: Calculate the Mass-Weighted Barysz Graph Energy for each molecule
    # and then identify the one with the lowest energy.
    molecule_energies = {}
    for name, molecule_symbols in y_molecules.items():
        n_atoms = len(molecule_symbols)
        if n_atoms <= 1:
            molecule_energies[name] = 0
            continue

        masses = [atomic_masses[sym] for sym in molecule_symbols]
        
        # Construct the Mass-weighted Laplacian matrix B, assuming a linear chain structure.
        B = np.zeros((n_atoms, n_atoms), dtype=float)
        
        # Diagonal elements are degree(i) / mass(i).
        B[0, 0] = 1.0 / masses[0]  # First atom, degree 1
        B[n_atoms - 1, n_atoms - 1] = 1.0 / masses[n_atoms - 1]  # Last atom, degree 1
        for i in range(1, n_atoms - 1):
            B[i, i] = 2.0 / masses[i]  # Middle atoms, degree 2
        
        # Off-diagonal elements for adjacent atoms are -1 / sqrt(mass(i) * mass(j)).
        for i in range(n_atoms - 1):
            val = -1.0 / np.sqrt(masses[i] * masses[i+1])
            B[i, i + 1] = val
            B[i + 1, i] = val
        
        # The energy is the sum of the absolute values of the eigenvalues of B.
        # np.linalg.eigvalsh is used as the matrix is symmetric.
        eigenvalues = np.linalg.eigvalsh(B)
        energy = np.sum(np.abs(eigenvalues))
        molecule_energies[name] = energy

    # Identify the molecule 'Y' with the minimum energy.
    identified_y_name = min(molecule_energies, key=molecule_energies.get)
    identified_y_energy = molecule_energies[identified_y_name]
    identified_y_symbols = y_molecules[identified_y_name]

    # Step 5 & 6: For the identified molecule Y, find its energy, and the min/max atomic masses.
    masses_of_identified_y = [atomic_masses[sym] for sym in identified_y_symbols]
    min_mass = min(masses_of_identified_y)
    max_mass = max(masses_of_identified_y)

    # Calculate the final requested product.
    result = identified_y_energy * min_mass * max_mass

    # Print the detailed results and the final answer in the required format.
    print(f"The element Y with the lowest Mass-Weighted Barysz Graph Energy is: {identified_y_name}")
    print(f"The interpreted molecular formula for {identified_y_name} is a chain: {'-'.join(identified_y_symbols)}")
    print(f"The calculated energy is: {identified_y_energy}")
    print(f"The minimum atomic mass in this molecule is: {min_mass}")
    print(f"The maximum atomic mass in this molecule is: {max_mass}")
    print("\nThe final result is the product of the energy, the minimum mass, and the maximum mass.")
    print(f"Final Calculation: {identified_y_energy} * {min_mass} * {max_mass} = {result}")
    print(f"<<<{result}>>>")

if __name__ == '__main__':
    solve_saber_duel_riddle()