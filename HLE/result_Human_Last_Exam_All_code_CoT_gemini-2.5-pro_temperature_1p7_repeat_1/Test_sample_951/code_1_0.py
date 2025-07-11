import sys
import io

# Save the original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

def predict_ir_phonons():
    """
    Predicts the number of IR-active phonons for LiNiPO4 in different polarizations.
    
    This function is based on the results of factor group analysis for crystals with
    the olivine structure, belonging to the Pnma space group (D2h point group).
    The analysis shows that for a unit cell with 28 atoms (Z=4 for LiNiPO4),
    the IR-active optical modes are distributed among the B1u, B2u, and B3u
    irreducible representations.
    """
    
    # Step 1: State the crystal properties
    crystal_system = "Orthorhombic"
    space_group = "Pnma (No. 62)"
    point_group = "D2h"
    formula_units_per_cell = 4
    atoms_per_cell = 28 # (1 Li + 1 Ni + 1 P + 4 O) * 4
    
    # Step 2: Define the number of IR-active optical modes from factor group analysis results
    # These are known results from literature for olivine structures (e.g., LiFePO4, LiNiPO4).
    # Total optical modes: Gamma_opt = 11Ag + 7B1g + 11B2g + 7B3g + 7Au + 10B1u + 6B2u + 10B3u
    ir_active_modes = {
        'B1u': 10,  # Number of B1u optical modes
        'B2u': 6,   # Number of B2u optical modes
        'B3u': 10,  # Number of B3u optical modes
    }
    
    # Step 3: Map polarizations (x,y,z) to crystallographic axes (a,b,c) and irreps
    # This is the standard convention for the Pnma setting.
    polarization_map = {
        'x': {'axis': 'a', 'irrep': 'B3u'},
        'y': {'axis': 'b', 'irrep': 'B2u'},
        'z': {'axis': 'c', 'irrep': 'B1u'},
    }
    
    # Step 4: Get the number of phonons for each polarization
    num_x = ir_active_modes[polarization_map['x']['irrep']]
    num_y = ir_active_modes[polarization_map['y']['irrep']]
    num_z = ir_active_modes[polarization_map['z']['irrep']]

    # Step 5: Print the result in the required format
    # The final equation shows the number of IR-active phonons for each polarization
    print(f"E||x: {num_x}, E||y: {num_y}, E||z: {num_z}")

predict_ir_phonons()

# Restore the original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the final result and the required answer tag
print(output.strip())
final_answer = output.strip()
print(f"<<<{final_answer}>>>")