import math

def calculate_exciton_binding_energy():
    """
    Calculates the binding energy of the n=3 exciton state in a 2D semiconductor.
    """
    # Given values from the problem description
    E_g = 3.0  # Band gap in eV
    E_1s_peak = 1.0  # 1s exciton resonance peak in eV
    n = 3      # Principal quantum number for the target state

    # Step 1: Calculate the 1s exciton binding energy (E_b1)
    # The energy of the exciton peak is E_peak = E_g - E_binding.
    # So, the binding energy is E_binding = E_g - E_peak.
    E_b1 = E_g - E_1s_peak

    # Step 2: Calculate the exciton Rydberg energy (R_X) for a 2D system
    # For a 2D system, the binding energy series is E_bn = R_X / (n - 0.5)^2.
    # For the ground state (n=1), this becomes E_b1 = R_X / (1 - 0.5)^2 = R_X / 0.25 = 4 * R_X.
    # Therefore, R_X = E_b1 / 4.
    R_X = E_b1 / 4.0

    # Step 3: Calculate the binding energy for the n=3 state (E_b3)
    # This is the "Rydberg energy for n = 3" requested by the user.
    denominator = (n - 0.5)**2
    E_b3 = R_X / denominator

    # Step 4: Print the results and the final equation with numerical values
    print("This script calculates the exciton binding energy for the n=3 state.\n")
    print(f"Given: Band Gap (E_g) = {E_g} eV, 1s Exciton Peak = {E_1s_peak} eV")
    print("-" * 50)
    
    print("Step 1: Calculate the 1s binding energy (E_b1)")
    print(f"E_b1 = E_g - E_1s_peak = {E_g} eV - {E_1s_peak} eV = {E_b1} eV\n")
    
    print("Step 2: Calculate the fundamental exciton Rydberg energy (R_X)")
    print(f"For a 2D system, E_b1 = R_X / (1 - 0.5)^2 = 4 * R_X")
    print(f"R_X = E_b1 / 4 = {E_b1} eV / 4 = {R_X} eV\n")

    print(f"Step 3: Calculate the binding energy for n={n} (E_b{n})")
    print("The final calculation is:")
    print(f"E_b{n} = R_X / (n - 0.5)^2")
    print(f"E_b{n} = {R_X} eV / ({n} - 0.5)^2")
    print(f"E_b{n} = {R_X} eV / ({n-0.5})^2")
    print(f"E_b{n} = {R_X} eV / {denominator}")
    print(f"E_b{n} = {E_b3:.2f} eV")

# Execute the function to print the solution
calculate_exciton_binding_energy()