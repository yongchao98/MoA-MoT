import sys

def find_optimal_ratio():
    """
    This function provides the experimentally-determined ideal Ni/Ce ratio
    for Ni-Ceria catalysts based on scientific literature.
    """
    
    # Based on numerous experimental studies, a well-cited optimal atomic ratio
    # that balances activity and stability is approximately 1 part Ni to 9 parts Ce.
    ni_parts = 1
    ce_parts = 9
    
    # Calculate the decimal representation of the ratio
    decimal_ratio = ni_parts / ce_parts
    
    # Print the explanation and results
    print("The ideal Ni/Ce ratio in Ni-Ceria nanoparticles is determined through experimental research, not a simple calculation.")
    print("The goal is to maximize the dispersion of active Ni sites and the strong metal-support interaction (SMSI) with the Ceria (CeO2) support.")
    print("\nBased on scientific literature, a frequently cited optimal atomic ratio is:")
    
    # Per the instruction to "output each number in the final equation!"
    print(f"Ni : Ce = {ni_parts} : {ce_parts}")
    
    print(f"\nThis corresponds to a decimal ratio of Ni/Ce = {decimal_ratio:.2f}.")
    print("\nThis composition ensures that Ni is highly dispersed as small nanoparticles or even single atoms on the ceria surface, which is crucial for maximizing catalytic activity and stability in the Water Gas Shift and water splitting reactions.")

if __name__ == "__main__":
    find_optimal_ratio()