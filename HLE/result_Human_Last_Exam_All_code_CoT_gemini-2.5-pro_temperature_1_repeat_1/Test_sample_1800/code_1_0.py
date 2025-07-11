def find_ideal_ratio():
    """
    This function explains and calculates a representative ideal Ni/Ce ratio
    for Ni-Ceria catalysts based on scientific literature.
    """
    # Based on research, a low Ni loading is optimal to maximize catalytic sites
    # and prevent particle agglomeration (sintering).
    # A frequently reported, highly effective molar ratio is 1 part Ni to 9 parts Ce.
    ni_parts = 1
    ce_parts = 9

    # Calculate the numerical ratio
    ideal_ratio = ni_parts / ce_parts

    print("The ideal Ni/Ce ratio in Ni-Ceria nanoparticles is a range rather than a single point, but a consensus in scientific literature points to a low nickel loading for optimal performance in WGS and WS reactions.")
    print("This maximizes the dispersion of Ni particles and the crucial interface with the Ceria support.")
    print("\nA commonly cited, highly effective molar ratio is:")
    print(f"{ni_parts} part Ni to {ce_parts} parts Ce")
    print("\nThe final equation for this ratio is:")
    # The final print statement includes each number from the calculation.
    print(f"{ni_parts} / {ce_parts} = {ideal_ratio:.2f}")

if __name__ == "__main__":
    find_ideal_ratio()