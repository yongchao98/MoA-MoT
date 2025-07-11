import math

def calculate_symmetry_factors():
    """
    Calculates and explains the symmetry factors for second-order vacuum bubble
    diagrams in phi^4 scalar field theory.
    """
    print("In phi^4 theory, there are two distinct connected second-order vacuum bubble diagrams.")
    print("We will calculate the symmetry factor for each and then sum them.")
    print("-" * 50)

    # --- Diagram 1: The "Melon" Diagram ---
    print("Diagram 1: The 'Melon' Diagram")
    print("Description: Two vertices connected directly by four identical propagators.")
    
    # Symmetries of the Melon diagram
    vert_perms = 2  # Swapping the two identical vertices
    prop_perms = math.factorial(4) # Permuting the four identical propagators
    
    s1 = vert_perms * prop_perms
    
    print("\nCalculation of symmetry factor S1:")
    print(f" - A factor of {vert_perms} from swapping the two identical vertices.")
    print(f" - A factor of {prop_perms} (which is 4!) from permuting the four identical propagators.")
    print(f"S1 = (vertex permutations) * (propagator permutations)")
    print(f"S1 = {vert_perms} * {prop_perms} = {s1}")
    print("-" * 50)

    # --- Diagram 2: The "Spectacles" Diagram ---
    print("Diagram 2: The 'Spectacles' Diagram")
    print("Description: Two vertices, each with one self-loop, connected by two propagators.")

    # Symmetries of the Spectacles diagram
    vert_swap = 2       # Swapping the two identical vertices
    prop_swap = 2       # Swapping the two propagators connecting the vertices
    loop1_flip = 2      # Flipping the direction of the first self-loop
    loop2_flip = 2      # Flipping the direction of the second self-loop

    s2 = vert_swap * prop_swap * loop1_flip * loop2_flip

    print("\nCalculation of symmetry factor S2:")
    print(f" - A factor of {vert_swap} from swapping the two identical vertices.")
    print(f" - A factor of {prop_swap} from swapping the two identical propagators connecting the vertices.")
    print(f" - A factor of {loop1_flip} from flipping the self-loop on the first vertex.")
    print(f" - A factor of {loop2_flip} from flipping the self-loop on the second vertex.")
    print(f"S2 = (vertex swap) * (propagator swap) * (loop 1 flip) * (loop 2 flip)")
    print(f"S2 = {vert_swap} * {prop_swap} * {loop1_flip} * {loop2_flip} = {s2}")
    print("-" * 50)
    
    # --- Summation ---
    total_sum = s1 + s2
    
    print("Sum of all symmetry factors:")
    # The final equation with each number explicitly printed
    print(f"{s1} + {s2} = {total_sum}")

if __name__ == '__main__':
    calculate_symmetry_factors()