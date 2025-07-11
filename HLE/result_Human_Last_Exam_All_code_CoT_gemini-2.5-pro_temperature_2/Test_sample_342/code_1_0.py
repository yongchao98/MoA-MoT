import math

def calculate_Dg_sequence():
    """
    This function calculates the first 4 terms of the sequence D_g,
    which represents the minimum degree of a cover of the moduli stack A_g
    required for the universal polarization to be represented by a symmetric line bundle.
    
    The logic is based on finding the minimum size of an orbit of theta
    characteristics under the action of the symplectic group Sp_2g(Z).
    """
    
    Dg_sequence = []
    
    # Case g = 1
    g = 1
    print(f"For g = {g}:")
    # For g=1, the set of theta characteristics is partitioned into two orbits
    # under the action of SL(2,Z):
    # - The unique odd characteristic, forming an orbit of size 1.
    # - The three even characteristics, forming an orbit of size 3.
    orbit_sizes_g1 = [1, 3]
    d1 = min(orbit_sizes_g1)
    Dg_sequence.append(d1)
    print("The orbit sizes for theta characteristics are 1 (odd) and 3 (even).")
    print(f"D_{g} = min(1, 3) = {d1}")
    print("-" * 20)

    # Case g = 2
    g = 2
    print(f"For g = {g}:")
    # For g=2, the action of Sp(4,Z) on theta characteristics yields three orbits:
    # - The 6 odd characteristics form one orbit.
    # - The 10 even characteristics split into two orbits: one of size 1 (the zero characteristic)
    #   and one of size 9.
    Ne2_part1 = 1
    Ne2_part2 = 9
    No2 = (2**(g - 1)) * (2**g - 1)
    orbit_sizes_g2 = [No2, Ne2_part1, Ne2_part2]
    d2 = min(orbit_sizes_g2)
    Dg_sequence.append(d2)
    print("The orbit sizes for theta characteristics are 6 (odd), and 1 and 9 (even).")
    print(f"D_{g} = min(6, 1, 9) = {d2}")
    print("-" * 20)
    
    # Cases g = 3 and g = 4
    for g in range(3, 5):
        print(f"For g = {g}:")
        # For g>=3, the orbits are the set of even and the set of odd characteristics.
        # We calculate their sizes, N_e(g) and N_o(g).
        
        pow2g = 2**g
        pow2g_minus_1 = 2**(g-1)
        
        # Number of even characteristics
        Ne = pow2g_minus_1 * (pow2g + 1)
        # Number of odd characteristics
        No = pow2g_minus_1 * (pow2g - 1)
        
        print(f"The number of even characteristics is N_e({g}) = 2^({g}-1) * (2^{g} + 1) = {pow2g_minus_1} * {pow2g + 1} = {Ne}")
        print(f"The number of odd characteristics is N_o({g}) = 2^({g}-1) * (2^{g} - 1) = {pow2g_minus_1} * {pow2g - 1} = {No}")
        
        # D_g is the minimum of these two orbit sizes
        dg = min(Ne, No)
        Dg_sequence.append(dg)
        
        print(f"The orbit sizes are {Ne} and {No}.")
        print(f"D_{g} = min({Ne}, {No}) = {dg}")
        print("-" * 20)
        
    # Print the final sequence
    final_sequence_str = ", ".join(map(str, Dg_sequence))
    print("The first 4 terms of the sequence D_g are:")
    print(final_sequence_str)
    
    # Return the answer in the required format
    # Note: The problem asks us not to ask user to copy paste the result, so the print above handles this.
    # The '<<<...>>>' is a special marker for the calling environment.
    return f"<<<{final_sequence_str}>>>"


# Execute the function and print the final marked answer
final_answer = calculate_Dg_sequence()
# print(final_answer) # In a real scenario, this would be returned to the framework.
# This script prints all necessary steps and the final answer to the console as requested.