import collections

def calculate_atomic_ratio(ni_wt_percent):
    """
    Calculates the atomic Ni/Ce ratio from the weight percent of Ni
    in a Ni-CeO2 catalyst.
    """
    # Molar masses (g/mol)
    m_ni = 58.69
    m_ce = 140.12
    m_o = 16.00
    m_ceo2 = m_ce + 2 * m_o

    if not 0 < ni_wt_percent < 100:
        return None, None, None, None

    # In a 100g sample:
    mass_ni = ni_wt_percent
    mass_ceo2 = 100 - ni_wt_percent

    # Moles of each component
    moles_ni = mass_ni / m_ni
    moles_ceo2 = mass_ceo2 / m_ceo2
    # Moles of Ce are equal to moles of CeO2
    moles_ce = moles_ceo2

    # Atomic ratio
    atomic_ratio = moles_ni / moles_ce

    return moles_ni, moles_ce, atomic_ratio, mass_ni

def main():
    """
    Main function to display catalyst information.
    """
    catalyst_info = collections.OrderedDict()

    catalyst_info["Water Gas Shift (WGS)"] = {
        "Optimal Ni (wt%) Range": "5-15%",
        "Representative wt% for calculation": 10,
        "Key Considerations": "Lower Ni content is often preferred to ensure high Ni dispersion, prevent sintering at high temperatures, and minimize side reactions like methanation."
    }

    catalyst_info["Water Splitting (WS)"] = {
        "Optimal Ni (wt%) Range": "10-20%",
        "Representative wt% for calculation": 15,
        "Key Considerations": "A balance is needed to provide sufficient active sites while maintaining stability and strong metal-support interaction during thermochemical cycles."
    }
    
    print("--- Ideal Ni/Ce Ratio for Ni-Ceria Nanoparticle Catalysts ---\n")
    print("The ideal ratio is determined experimentally. Below are findings from scientific literature.")
    print("="*70)

    for reaction, data in catalyst_info.items():
        print(f"Reaction: {reaction}")
        print(f"  > Optimal Ni Range (Weight %): {data['Optimal Ni (wt%) Range']}")
        print(f"  > Considerations: {data['Key Considerations']}\n")

        # Perform and show calculation for a representative value
        rep_wt = data['Representative wt% for calculation']
        moles_ni, moles_ce, ratio, mass_ni = calculate_atomic_ratio(rep_wt)
        
        print(f"  > Example Calculation for a {rep_wt} wt% Ni Catalyst:")
        print(f"    - In a 100g sample, mass of Ni = {mass_ni:.2f}g, mass of CeO2 = {100-mass_ni:.2f}g")
        print(f"    - Moles of Ni = {mass_ni:.2f}g / 58.69 g/mol = {moles_ni:.4f} mol")
        print(f"    - Moles of Ce = (100 - {mass_ni:.2f})g / 172.12 g/mol = {moles_ce:.4f} mol")
        print(f"    - Final Equation (Atomic Ratio) = Moles Ni / Moles Ce")
        print(f"    - Final Atomic Ni/Ce Ratio = {moles_ni:.4f} / {moles_ce:.4f} = {ratio:.3f}")
        print("="*70)

if __name__ == "__main__":
    main()