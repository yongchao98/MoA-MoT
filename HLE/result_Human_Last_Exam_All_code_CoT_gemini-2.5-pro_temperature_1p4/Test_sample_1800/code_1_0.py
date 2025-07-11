import sys

def find_ideal_ratio():
    """
    This function explains and provides the ideal Ni/Ce atomic ratio
    for catalytic applications like WGS and water splitting based on
    a review of scientific literature.
    """
    # Explanation based on scientific consensus
    print("The ideal Ni/Ce ratio in Ni-Ceria catalysts is not a single universal value,")
    print("as it depends on the synthesis method and specific reaction conditions.")
    print("However, extensive research has shown that a high catalytic performance is typically")
    print("achieved with a relatively low nickel loading. This ensures high dispersion of Ni")
    print("nanoparticles and maximizes the crucial interaction with the Ceria support,")
    print("which enhances the creation of oxygen vacancies vital for the reactions.\n")

    # Define the parts of the commonly cited optimal atomic ratio
    ni_part = 1
    ce_part = 9

    # The final "equation" as requested
    print("Based on numerous studies, a highly effective atomic ratio of Ni to Ce is:")
    print(f"    Ni : Ce = {ni_part} : {ce_part}")

    # Corresponding molar ratio value
    ratio_value = ni_part / ce_part
    print(f"\nThis corresponds to a Ni/Ce molar ratio of approximately {ratio_value:.2f}.")

if __name__ == '__main__':
    find_ideal_ratio()