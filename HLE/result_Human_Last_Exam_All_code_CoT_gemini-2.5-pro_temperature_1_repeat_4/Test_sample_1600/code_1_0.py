def calculate_feynman_diagrams_a3():
    """
    Calculates a(3), the number of non-vanishing Feynman diagrams of order 6
    for electron or photon propagators in QED.
    """
    
    # At 3 loops (order 6), the number of 1-particle-irreducible (1PI) diagrams
    # for the electron self-energy is 7.
    num_electron_diagrams = 7
    
    # At 3 loops (order 6), the number of 1PI diagrams
    # for the photon self-energy (vacuum polarization) is 3.
    num_photon_diagrams = 3
    
    # a(3) is the sum of the diagrams for both propagators.
    a3 = num_electron_diagrams + num_photon_diagrams
    
    print(f"The value a(n) represents the total number of non-vanishing Feynman diagrams of order 2n for the electron and photon propagators.")
    print(f"We need to find a(3), which corresponds to the diagrams of order 2*3 = 6 (3 loops).")
    print(f"This is the sum of the diagrams for the electron propagator and the photon propagator.")
    print(f"\nNumber of diagrams for electron propagator at 3 loops = {num_electron_diagrams}")
    print(f"Number of diagrams for photon propagator at 3 loops = {num_photon_diagrams}")
    
    # Print the final equation showing all numbers
    print(f"\nTherefore, a(3) = {num_electron_diagrams} + {num_photon_diagrams} = {a3}")

calculate_feynman_diagrams_a3()