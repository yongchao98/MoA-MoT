import numpy as np

def demonstrate_polarization_invariance():
    """
    This function demonstrates that a system of scalar operators cannot
    generate a vector beam from a linearly polarized input.
    """
    # Let's represent the beam cross-section as a 2D grid
    grid_size = 100

    # 1. Define the Input Beam
    # The input beam has a complex spatial profile but is uniformly linearly polarized.
    # We model this by giving it a non-zero Ex component and a zero Ey component.
    x = np.linspace(-1, 1, grid_size)
    y = np.linspace(-1, 1, grid_size)
    X, Y = np.meshgrid(x, y)
    
    # A simple Gaussian spatial profile for the Ex component
    ex_component_in = np.exp(-(X**2 + Y**2))
    # The Ey component is zero everywhere
    ey_component_in = np.zeros((grid_size, grid_size), dtype=complex)

    print("--- Analysis of the Optical System ---")
    print(f"Initial state: Does the orthogonal (Ey) component exist? {np.any(ey_component_in != 0)}")

    # 2. Define the System's Scalar Operators
    # These operators act on the spatial profile but not the polarization.
    # We define them to act on a single component (like Ex or Ey).

    # Operator for a phase-shaping medium (e.g., a random phase mask)
    random_phase_mask = np.exp(1j * 2 * np.pi * np.random.rand(grid_size, grid_size))
    def apply_phase_mask(field_component):
        return field_component * random_phase_mask

    # Operator for free-space propagation (simplified as a Fourier transform)
    def propagate_in_free_space(field_component):
        # The key is that this operation is scalar; it acts on a single component.
        return np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(field_component)))

    # 3. Propagate the beam through the system
    # We apply the operators sequentially to each polarization component.
    
    # Current state of the components
    ex_current = ex_component_in
    ey_current = ey_component_in

    # Path 1: Free space -> Phase mask -> Free space
    ex_current = propagate_in_free_space(ex_current)
    ey_current = propagate_in_free_space(ey_current) # Propagating zero results in zero

    ex_current = apply_phase_mask(ex_current)
    ey_current = apply_phase_mask(ey_current) # Multiplying zero by a mask results in zero

    ex_current = propagate_in_free_space(ex_current)
    ey_current = propagate_in_free_space(ey_current) # Propagating zero results in zero
    
    # The problem describes a second path, but the logic remains the same.
    # Any further scalar operation (like passing through an inverse mask) will also
    # keep the Ey component at zero if it's already zero.

    # 4. Analyze the Final Output
    ex_final = ex_current
    ey_final = ey_current

    print(f"Final state: Does the orthogonal (Ey) component exist? {np.any(ey_final != 0)}")
    
    # 5. Conclusion
    print("\n--- Conclusion ---")
    print("A vector beam requires a spatially varying polarization, which means it must have non-zero orthogonal components (e.g., Ex and Ey) with different spatial profiles.")
    print("The simulation shows that if the input beam is purely linearly polarized (Ey = 0), and the system only contains scalar operators (like phase masks and free-space propagation), the Ey component remains zero throughout.")
    print("Therefore, the system described cannot generate a vector beam from a linearly polarized input.")
    print("\nThe answer to the question 'Can we get an arbitrary vector beam...?' is:")
    print("No")

demonstrate_polarization_invariance()