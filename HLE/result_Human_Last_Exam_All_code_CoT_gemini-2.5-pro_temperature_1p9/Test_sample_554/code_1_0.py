import sys

def calculate_fluorine_atoms():
    """
    Calculates the number of fluorine atoms in a hypothetical perfluoronanocar
    based on the structure of the first 4-wheeled nanocar.
    """
    # Plan:
    # 1. A "perfluoronanocar" is a nanocar where all hydrogen (H) atoms on the
    #    hydrocarbon chassis are replaced by fluorine (F) atoms.
    # 2. We use the structure of the original 4-wheeled nanocar (Tour group, 2005).
    # 3. The hydrocarbon framework consists of a central benzene ring and four "arms".
    # 4. We count the hydrogens on this framework to find the number of fluorines.

    # Number of hydrogens on the central benzene ring of the chassis.
    # The ring is 1,2,4,5-tetrasubstituted (C6H2), so it has 2 hydrogens.
    hydrogens_central_ring = 2

    # Number of arms connecting the chassis to the wheels.
    num_arms = 4
    
    # Each arm contains one para-substituted phenylene group (C6H4), which has 4 hydrogens.
    hydrogens_per_arm_phenylene = 4

    # The "perfluoro-" prefix means the total number of fluorine atoms will equal
    # the total number of hydrogen atoms in the original structure.
    
    # Print out the step-by-step calculation.
    print("To find the number of fluorine atoms, we first count the replaceable hydrogen atoms in a standard nanocar's chassis:")
    print(f"1. Hydrogens on the central benzene ring: {hydrogens_central_ring}")
    print(f"2. Hydrogens on the phenylene groups in the four arms: {num_arms} arms * {hydrogens_per_arm_phenylene} hydrogens per arm")
    
    print("\nThe final equation sums these parts:")
    print(f"Total Fluorine Atoms = (Hydrogens on Central Ring) + (Number of Arms) * (Hydrogens per Arm)")
    
    # Explicitly showing the numbers in the equation
    print(f"Total Fluorine Atoms = {hydrogens_central_ring} + {num_arms} * {hydrogens_per_arm_phenylene}")
    
    # Perform the calculation
    hydrogens_on_arms = num_arms * hydrogens_per_arm_phenylene
    total_fluorines = hydrogens_central_ring + hydrogens_on_arms
    
    print(f"Total Fluorine Atoms = {hydrogens_central_ring} + {hydrogens_on_arms} = {total_fluorines}")
    
    print(f"\nTherefore, a hypothetical perfluoronanocar would contain {total_fluorines} fluorine atoms.")

calculate_fluorine_atoms()