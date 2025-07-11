import math

def solve_exciton_energy():
    """
    Calculates the Rydberg energy for the n=3 exciton state in a 2D semiconductor.
    """
    # Given parameters
    Eg = 3.0  # Band gap in eV
    E_1s = 1.0 # 1s exciton resonance peak in eV
    n = 3     # Principal quantum number for the target state

    print("--- Step 1: Calculate the 1s exciton binding energy (Eb_1s) ---")
    # The exciton energy level E_n is related to the band gap Eg and binding energy Eb_n
    # by the formula: E_n = Eg - Eb_n.
    # Therefore, Eb_n = Eg - E_n.
    Eb_1s = Eg - E_1s
    print(f"The 1s exciton binding energy is calculated as:")
    print(f"Eb_1s = Eg - E_1s = {Eg} eV - {E_1s} eV = {Eb_1s} eV\n")


    print("--- Step 2: Calculate the effective Rydberg energy constant (R_X) ---")
    # The binding energy series for an ideal 2D exciton is Eb_n = R_X / (n - 1/2)^2.
    # For the ground state (n=1), this becomes Eb_1s = R_X / (1 - 0.5)^2 = R_X / (0.5)^2 = 4 * R_X.
    # From this, we can find R_X.
    R_X = Eb_1s / 4.0
    print("For a 2D system, Eb_n = R_X / (n - 0.5)^2. For n=1, Eb_1s = 4 * R_X.")
    print(f"The effective Rydberg constant is calculated as:")
    print(f"R_X = Eb_1s / 4 = {Eb_1s} eV / 4.0 = {R_X} eV\n")


    print("--- Step 3: Calculate the Rydberg energy for n=3 (Eb_3) ---")
    # Now we use the formula for n=3 to find its binding energy.
    # The question "Rydberg energy for n=3" refers to this value.
    denominator_val = (n - 0.5)
    denominator_sq = denominator_val**2
    Eb_3 = R_X / denominator_sq
    print(f"The Rydberg energy for n={n} is calculated using Eb_{n} = R_X / (n - 0.5)^2:")
    print(f"Eb_{n} = {R_X} eV / ({n} - 0.5)^2 = {R_X} eV / ({denominator_val})^2 = {R_X} eV / {denominator_sq} = {Eb_3:.3f} eV\n")
    
    print(f"The final answer for the Rydberg energy for n = 3 is {Eb_3:.3f} eV.")
    
    # Final answer in the required format
    print(f"\n<<<{Eb_3:.3f}>>>")

solve_exciton_energy()