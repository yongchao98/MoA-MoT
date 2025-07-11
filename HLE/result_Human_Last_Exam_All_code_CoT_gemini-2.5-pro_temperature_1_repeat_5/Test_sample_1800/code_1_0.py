import sys

def find_ideal_ratio():
    """
    This function presents the ideal Ni/Ce atomic ratio for catalytic applications
    based on findings from scientific literature.
    """
    # Based on research, an optimal atomic ratio is often cited where Nickel
    # constitutes about 10% of the metal atoms. This is represented by a 1:9 ratio.
    # This composition maximizes the perimeter interface between Ni and Ceria,
    # which is crucial for the catalytic activity in WGS and WS reactions.
    ni_parts = 1
    ce_parts = 9

    # Calculate the total atomic parts
    total_parts = ni_parts + ce_parts

    # Calculate the atomic percentage of Nickel
    ni_atomic_percent = (ni_parts / total_parts) * 100

    print("The ideal Ni/Ce ratio is one that maximizes the active Ni-Ceria interface.")
    print("Based on numerous studies, a highly effective atomic ratio is found to be:")
    print(f"Ni : Ce = {ni_parts} : {ce_parts}")
    print("\nThis ratio corresponds to a specific atomic loading of Nickel.")
    print("The calculation for the atomic percentage of Ni is:")
    print(f"Ni % = (Ni_parts / (Ni_parts + Ce_parts)) * 100")
    print(f"Ni % = ({ni_parts} / ({ni_parts} + {ce_parts})) * 100 = {ni_atomic_percent:.1f}%")
    print("\nTherefore, an ideal catalyst composition is approximately 10% atomic Nickel relative to Cerium.")

if __name__ == '__main__':
    find_ideal_ratio()
