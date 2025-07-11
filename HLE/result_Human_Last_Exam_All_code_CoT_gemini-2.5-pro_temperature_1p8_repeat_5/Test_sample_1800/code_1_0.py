import sys

def find_ideal_ratio():
    """
    Calculates and explains the ideal Ni/Ce molar ratio for catalysis
    based on findings from scientific literature.
    """
    
    # --- Explanation ---
    print("Finding the Ideal Ni/Ce Ratio for Catalysis")
    print("-" * 50)
    print("For Ni-Ceria catalysts used in the Water Gas Shift (WGS) and Water Splitting (WS) reactions, the performance is highly dependent on the ratio of Nickel (Ni) to Cerium (Ce).")
    print("\nScientific Background:")
    print("1. Ceria (CeO₂) is an excellent support material due to its high oxygen storage capacity (OSC).")
    print("2. Nickel (Ni) is the active metal that catalyzes the reaction.")
    print("3. The optimal ratio creates a strong synergy, forming a Ni-Ce-O solid solution at the interface. This enhances catalytic activity by keeping Ni particles small and stable, maximizing the active sites.")
    print("   - Too little Ni results in insufficient active sites.")
    print("   - Too much Ni leads to the formation of larger NiO clusters, which are less active and prone to deactivation.")

    # --- Calculation based on Scientific Literature ---
    print("\nCalculation from Literature Data:")
    print("Research frequently indicates that the optimal performance is achieved when the catalyst forms a solid solution with the formula NiₓCe₁-ₓO₂.")
    
    # A value for x around 0.2 is commonly cited as being optimal for creating
    # a highly active and stable solid solution.
    optimal_x = 0.2
    
    print(f"A commonly reported optimal value for 'x' in the formula NiₓCe₁-ₓO₂ is {optimal_x}.")
    
    # In one mole of the NiₓCe₁-ₓO₂ compound:
    moles_ni = optimal_x
    moles_ce = 1 - optimal_x
    
    print(f"\nThis corresponds to:")
    print(f" - Molar part of Ni = {moles_ni}")
    print(f" - Molar part of Ce = 1 - {optimal_x} = {moles_ce:.1f}")

    # Calculate the ideal ratio
    ideal_ratio = moles_ni / moles_ce
    
    print("\nThe ideal Ni/Ce molar ratio is the molar part of Ni divided by the molar part of Ce.")
    print("\nFinal Equation and Result:")
    # Print out each number in the final equation as requested
    print(f"The ideal Ni/Ce molar ratio is {moles_ni} / {moles_ce:.1f} = {ideal_ratio}")
    print("\nThis means the ideal catalyst has a molar ratio of 1 part Nickel to 4 parts Cerium (1:4).")
    print("-" * 50)
    
    # Storing the final numerical answer for the platform.
    # Note: Using sys.stdout for clarity on where the final 'answer' goes.
    sys.stdout.write(f"\n<<<ANSER_BLOCK_FORMAT>>>\n")
    sys.stdout.write(f"<<<{ideal_ratio}>>>")


if __name__ == '__main__':
    find_ideal_ratio()
