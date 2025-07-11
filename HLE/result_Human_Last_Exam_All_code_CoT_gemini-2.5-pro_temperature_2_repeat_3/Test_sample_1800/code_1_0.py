def analyze_ni_ceria_ratio():
    """
    Analyzes and presents the ideal Ni/Ce ratio for catalytic reactions
    based on findings from scientific literature.
    """
    print("Analysis of Ideal Ni/Ce Ratio in Ni-Ceria Catalysts")
    print("-" * 55)

    print("\n[Background]")
    print("The goal is to maximize the synergy between Nickel (the active catalyst) and Ceria (the support).")
    print("Too little Ni leads to insufficient active sites.")
    print("Too much Ni can lead to particle agglomeration (sintering), reducing efficiency and stability.")

    print("\n[1. Water Gas Shift (WGS) Reaction]")
    print("For WGS, a high dispersion of Ni is critical to prevent side reactions like methanation.")
    print("Literature often points to an optimal range where Ni content is significant but not excessive.")
    print("An effective Ni to Ce molar ratio is often found in the range of 1:9 to 1:4.")
    print("\nLet's break down the numbers in this common optimal range:")
    print("Lower End of the Range (Ratio 1:9):")
    print("  - Nickel (Ni) parts: 1")
    print("  - Cerium (Ce) parts: 9")
    print("  (This corresponds to approximately 10 mol% Nickel)")
    
    print("\nHigher End of the Range (Ratio 1:4):")
    print("  - Nickel (Ni) parts: 1")
    print("  - Cerium (Ce) parts: 4")
    print("  (This corresponds to approximately 20 mol% Nickel)")


    print("\n[2. Water Splitting (WS)]")
    print("Similar to WGS, the interface between Ni and Ceria is crucial for activating water molecules.")
    print("The optimal ratio for water splitting often falls within a similar range as WGS.")
    print("High performance is consistently reported for catalysts with Ni content between 10-20 mol%.")
    print("This again corresponds to a Ni:Ce molar ratio of approximately 1:9 to 1:4, which ensures well-dispersed NiO nanoparticles on the CeO2 support.")
    
    print("\n[Conclusion]")
    print("Across both WGS and WS applications, a molar ratio of Ni to Ce from 1:9 to 1:4 is a widely cited ideal range to maximize catalytic activity and stability.")

if __name__ == '__main__':
    analyze_ni_ceria_ratio()