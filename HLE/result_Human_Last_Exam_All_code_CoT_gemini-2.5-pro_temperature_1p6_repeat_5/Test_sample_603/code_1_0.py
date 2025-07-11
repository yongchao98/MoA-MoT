import math

def calculate_dna_tm(sequence, oligo_conc_M=2e-7):
    """
    Calculates the melting temperature (Tm) of a DNA duplex using the
    Nearest-Neighbor Base-Pair (NNBP) model.

    This code uses thermodynamic parameters from SantaLucia, J. (1998), PNAS.
    """

    # Thermodynamic parameters for DNA nearest-neighbors in 1M NaCl
    # dH in kcal/mol, dS in cal/mol-K
    dH_table = {
        'AA': -7.6, 'TT': -7.6, 'AT': -7.2, 'TA': -7.2, 'CA': -8.5,
        'TG': -8.5, 'GT': -8.4, 'AC': -8.4, 'CT': -7.8, 'AG': -7.8,
        'GA': -8.2, 'TC': -8.2, 'CG': -10.6, 'GC': -9.8, 'GG': -8.0, 'CC': -8.0
    }

    dS_table = {
        'AA': -21.3, 'TT': -21.3, 'AT': -20.4, 'TA': -21.3, 'CA': -22.7,
        'TG': -22.7, 'GT': -22.4, 'AC': -22.4, 'CT': -21.0, 'AG': -21.0,
        'GA': -22.2, 'TC': -22.2, 'CG': -27.2, 'GC': -24.4, 'GG': -19.9, 'CC': -19.9
    }
    
    # Helix initiation corrections
    dH_init = 0.0  # kcal/mol
    dS_init = -1.4 # cal/mol-K (for non-self-complementary sequences)

    total_dH = dH_init
    total_dS = dS_init

    # Summing the thermodynamic values for each nearest-neighbor pair
    for i in range(len(sequence) - 1):
        pair = sequence[i:i+2]
        total_dH += dH_table.get(pair, 0)
        total_dS += dS_table.get(pair, 0)
        
    # Ideal gas constant in cal/mol-K
    R = 1.987 

    # Melting Temperature (Tm) formula for non-self-complementary duplexes:
    # Tm (°C) = (Total dH * 1000) / (Total dS + R * ln(Oligo Concentration)) - 273.15
    # The concentration term accounts for the entropy of bringing two strands together.
    tm_celsius = (total_dH * 1000) / (total_dS + R * math.log(oligo_conc_M)) - 273.15

    print(f"Analysis for DNA sequence: {sequence}\n")
    print("--- Equation for Melting Temperature (Tm) ---")
    print("Tm(°C) = (ΔH° * 1000) / (ΔS° + R * ln(C)) - 273.15\n")
    print("--- Values for the final equation ---")
    print(f"Total Enthalpy (ΔH°): {total_dH:.2f} kcal/mol")
    print(f"Total Entropy (ΔS°): {total_dS:.2f} cal/mol-K")
    print(f"Ideal Gas Constant (R): {R} cal/mol-K")
    print(f"Oligo Concentration (C): {oligo_conc_M} M")
    print("------------------------------------------")
    print(f"\nCalculated Melting Temperature (Tm): {tm_celsius:.2f} °C")


if __name__ == '__main__':
    # Example non-self-complementary DNA sequence
    dna_sequence = "GGTATACC"
    calculate_dna_tm(dna_sequence)