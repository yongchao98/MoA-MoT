import sys

def calculate_ideal_ratio():
    """
    Calculates and explains the ideal Ni/Ce ratio for Ni-Ceria catalysts
    based on values commonly reported in scientific literature.
    """
    # The ideal Ni/Ce ratio is determined experimentally. A review of literature
    # suggests that a nickel concentration of around 10 atomic percent is often
    # optimal for maximizing catalytic performance in Water Gas Shift (WGS) and
    # water splitting reactions.
    # This means for every 100 total metal atoms, 10 are Ni and 90 are Ce.

    ni_parts = 10
    ce_parts = 90

    # The Ni/Ce molar ratio is calculated by dividing the moles (or atomic parts)
    # of Ni by the moles (or atomic parts) of Ce.
    ideal_ratio = ni_parts / ce_parts

    print("Step 1: Identify the optimal atomic composition from scientific literature.")
    print(f"A common optimal composition is {ni_parts} atomic percent Ni and {ce_parts} atomic percent Ce.\n")

    print("Step 2: Express this composition as a molar ratio of Ni to Ce.")
    print("The equation is: Ratio = (Parts of Ni) / (Parts of Ce)\n")

    print("Step 3: Calculate the final ratio using the values.")
    # We use sys.stdout.write to print the equation without a newline at the end
    # to fulfill the requirement of showing each number in the final equation.
    sys.stdout.write(f"Final Calculation: {ni_parts} / {ce_parts} = ")
    print(f"{ideal_ratio:.3f}\n") # Using print here adds the newline after the result

    print("Conclusion:")
    print(f"An ideal Ni/Ce molar ratio for Ni-Ceria catalysts is approximately {ideal_ratio:.3f}.")

if __name__ == '__main__':
    calculate_ideal_ratio()
