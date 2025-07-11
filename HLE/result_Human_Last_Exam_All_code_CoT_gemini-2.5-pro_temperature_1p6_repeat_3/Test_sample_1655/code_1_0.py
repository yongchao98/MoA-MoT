import numpy as np

def analyze_beam_transformation():
    """
    This script conceptually demonstrates why a scalar beam cannot be converted
    into a vector beam using only isotropic optical elements like phase masks
    and free space propagation.
    """

    # --- Setup ---
    grid_size = 128

    # 1. Define the Input Beam
    # We create a simple Gaussian beam profile as a stand-in for the Laguerre-Gaussian mode.
    # The key is that it is linearly polarized, so its Ey component is zero.
    x = np.linspace(-3, 3, grid_size)
    X, Y = np.meshgrid(x, y)
    spatial_profile = np.exp(-(X**2 + Y**2))
    
    # The beam is represented by its polarization components.
    # For a linearly, horizontally polarized beam, Ey is zero everywhere.
    input_beam = {
        'Ex': spatial_profile.astype(complex),
        'Ey': np.zeros_like(spatial_profile, dtype=complex)
    }

    print("--- Analysis of Optical System ---")
    print("\nStep 1: Define the Input Beam")
    print("The input beam has fixed linear polarization (it's a 'scalar beam').")
    print(f"Does the input beam have a non-zero Ey component? {not np.all(input_beam['Ey'] == 0)}")


    # 2. Define the Optical Elements (as functions)

    def free_space_prop(field_component):
        """Dummy function for free space propagation (e.g., via Fourier Transform).
           It acts on each polarization component independently."""
        # A real propagator is complex, but for this demo, the key is that the
        # function is called independently on Ex and Ey. We'll just return the input.
        return field_component

    def apply_phase_mask(field_component, mask):
        """Dummy function for a phase-shaping medium.
           It acts on each polarization component independently."""
        return field_component * mask

    # Create a random phase mask 'T' and its inverse.
    random_phase = np.random.rand(grid_size, grid_size) * 2 * np.pi
    T_mask = np.exp(1j * random_phase)
    T_inverse_mask = np.conj(T_mask)

    # 3. Simulate the full beam path

    # First system: Free Space -> T -> Free Space
    propagated_ex_1 = free_space_prop(input_beam['Ex'])
    propagated_ey_1 = free_space_prop(input_beam['Ey']) # This remains zero

    masked_ex_1 = apply_phase_mask(propagated_ex_1, T_mask)
    masked_ey_1 = apply_phase_mask(propagated_ey_1, T_mask) # This remains zero

    output1_ex = free_space_prop(masked_ex_1)
    output1_ey = free_space_prop(masked_ey_1) # This remains zero

    print("\nStep 2: Trace Beam to 'Output 1'")
    print(f"After passing through the first system, does the beam have a non-zero Ey component? {not np.all(output1_ey == 0)}")

    # Second system: Output 1 -> Free Space -> T_inverse
    propagated_ex_2 = free_space_prop(output1_ex)
    propagated_ey_2 = free_space_prop(output1_ey) # This remains zero

    output2_ex = apply_phase_mask(propagated_ex_2, T_inverse_mask)
    output2_ey = apply_phase_mask(propagated_ey_2, T_inverse_mask) # This remains zero
    
    print("\nStep 3: Trace Beam to 'Output 2'")
    print(f"After passing through the second system, does the final beam have a non-zero Ey component? {not np.all(output2_ey == 0)}")
    
    # 4. Conclusion
    print("\n--- Conclusion ---")
    print("A vector beam requires a spatially varying polarization, which means both its Ex and Ey components must be non-zero.")
    print("As demonstrated, since the input beam had Ey=0 and all optical elements were isotropic (acting on Ex and Ey independently), the Ey component remained zero throughout.")
    print("Therefore, the system cannot generate a vector beam from the given input.")
    
    print("\nFinal Answer to the Question:")
    # The final answer to the conceptual question.
    print("No")


# Run the analysis function
analyze_beam_transformation()