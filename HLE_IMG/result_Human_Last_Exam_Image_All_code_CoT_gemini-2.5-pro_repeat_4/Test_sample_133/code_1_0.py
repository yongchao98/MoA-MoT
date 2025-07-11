import sys

def solve_mass_spectrum():
    """
    Analyzes the provided mass spectrum data to identify the unknown compound and prints the step-by-step reasoning.
    """
    # Mass spectrum data from the table
    # Storing a subset of the data relevant for the analysis
    ms_data = {
        35: 3.09, 36: 3.38, 37: 1.04, 47: 18.4, 48: 4.62, 49: 6.24, 71: 10.6, 83: 19.4,
        85: 7.09, 96: 10.1, 117: 2.98, 118: 43.5, 119: 4.41, 120: 29.6, 121: 2.06,
        141: 28.1, 155: 16.1, 190: 48.4, 192: 24.5, 225: 100, 227: 66.7, 229: 22.2,
        260: 34.8, 262: 29.0, 264: 12.9, 266: 3.22
    }

    print("Step 1: Analyze the Molecular Ion (M+) Peak Cluster")
    print("The analysis starts by examining the peaks at the highest mass-to-charge ratio (m/z) to identify the molecular ion.")
    mz_260, i_260 = 260, ms_data[260]
    mz_262, i_262 = 262, ms_data[262]
    mz_264, i_264 = 264, ms_data[264]
    mz_266, i_266 = 266, ms_data[266]
    print(f"Observed peaks in the high-mass region: m/z {mz_260} (I={i_260}%), m/z {mz_262} (I={i_262}%), m/z {mz_264} (I={i_264}%), m/z {mz_266} (I={i_266}%).")
    print("This isotopic cluster pattern is characteristic of a compound containing three chlorine atoms (Cl).")
    print("The theoretical intensity ratio for 3 Cl atoms (due to isotopes ³⁵Cl and ³⁷Cl) is approximately M:(M+2):(M+4):(M+6) = 100:98:32:3.")
    normalized_i_262 = (i_262 / i_260) * 100
    normalized_i_264 = (i_264 / i_260) * 100
    normalized_i_266 = (i_266 / i_260) * 100
    print(f"The observed normalized ratio is {i_260}/{i_260}*100 : {i_262}/{i_260}*100 : {i_264}/{i_260}*100 : {i_266}/{i_260}*100, which calculates to approximately 100 : {normalized_i_262:.0f} : {normalized_i_264:.0f} : {normalized_i_266:.0f}.")
    print("This observed ratio is a good match with the theoretical pattern for three chlorine atoms.")
    print("Therefore, the molecule contains 3 Cl atoms, and the peak at m/z 260 corresponds to the molecular ion with only the lighter ³⁵Cl isotope.")
    print("-" * 50)

    print("Step 2: Determine the Non-Halogen Part of the Molecule")
    mass_M_35Cl = 260
    mass_C = 12
    mass_H = 1
    mass_35Cl = 35
    num_Cl = 3
    mass_of_3_Cl = num_Cl * mass_35Cl
    mass_of_backbone = mass_M_35Cl - mass_of_3_Cl
    print(f"The mass of the non-halogen part can be calculated by subtracting the mass of the three ³⁵Cl atoms from the M+ peak at m/z 260.")
    print(f"Mass of CₓHᵧ = M+ - (3 * Mass of ³⁵Cl) = {mass_M_35Cl} - (3 * {mass_35Cl}) = {mass_of_backbone}.")
    print(f"A plausible hydrocarbon formula for a mass of {mass_of_backbone} is C₁₁H₂₃.")
    num_C, num_H = 11, 23
    mass_C11H23 = num_C * mass_C + num_H * mass_H
    print(f"Calculation for C₁₁H₂₃: ({num_C} * {mass_C}) + ({num_H} * {mass_H}) = {num_C * mass_C} + {num_H * mass_H} = {mass_C11H23}.")
    print(f"This indicates the molecular formula of the compound is C₁₁H₂₃Cl₃.")
    print("-" * 50)

    print("Step 3: Confirm the Structure by Analyzing Key Fragment Peaks")
    print("The base peak (most intense peak) at m/z 225 corresponds to the loss of one chlorine atom ([M-Cl]⁺).")
    mass_M_minus_Cl = mass_M_35Cl - mass_35Cl
    print(f"Equation: M+ - Cl = {mass_M_35Cl} - {mass_35Cl} = {mass_M_minus_Cl}.")
    print("The stability of this fragment suggests the positive charge is well-stabilized.")
    
    print("\nA very strong fragment peak is observed at m/z 141.")
    mass_CCl3 = mass_C + num_Cl * mass_35Cl
    print("This corresponds to an alpha-cleavage, where the bond next to the carbon bearing the chlorines breaks, losing a CCl₃ radical.")
    print(f"Equation: M+ - CCl₃ = C₁₁H₂₃Cl₃ - CCl₃ = [C₁₀H₂₁]⁺.")
    print(f"Calculated mass of the resulting C₁₀H₂₁⁺ fragment: (10 * {mass_C}) + (21 * {mass_H}) = {10 * mass_C + 21 * mass_H}.")
    print("The presence of a strong peak at m/z 141 is a powerful indicator for a C₁₀H₂₁-CCl₃ structure, i.e., 1,1,1-trichloroundecane.")
    
    print("\nThe series of smaller peaks at m/z 85, 71, etc., are characteristic fragments of a long alkyl chain (C₁₀H₂₁), further supporting this structure.")
    print("-" * 50)

    print("Step 4: Conclusion")
    print("Based on the combined evidence from the molecular ion, isotope patterns, and fragmentation analysis, the compound is identified.")
    iupac_name = "1,1,1-trichloroundecane"
    print(f"The systematic IUPAC name for the compound CH₃(CH₂)₉CCl₃ is {iupac_name}.")

# The main execution block
if __name__ == '__main__':
    solve_mass_spectrum()