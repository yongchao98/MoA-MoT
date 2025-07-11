import numpy as np

def demonstrate_scalar_operation():
    """
    This function demonstrates that a scalar optical system cannot generate a
    vector beam from a uniformly polarized input beam.
    """
    # A vector beam is a beam with a spatially non-uniform polarization state.
    # To create one from a uniformly polarized beam, the system must be able
    # to manipulate polarization, i.e., it must be a non-scalar system
    # (like a q-plate or a birefringent medium).
    
    # The optical elements in the question are free space and a "phase-shaping medium".
    # Both of these are scalar operators, meaning they act on the phase and
    # amplitude of the beam but do not mix polarization components.

    print("Step 1: Define the input beam.")
    # Let's represent the input beam's electric field on a small 2x2 grid for demonstration.
    # The beam has two orthogonal polarization components, Ex and Ey.
    # The input beam is uniformly and linearly polarized, so we set the Ey component to zero.
    Ex_in = np.array([[1.0, 1.2], [1.3, 0.9]])
    Ey_in = np.array([[0.0, 0.0], [0.0, 0.0]])
    print("Input Ex component:\n", Ex_in)
    print("Input Ey component:\n", Ey_in)
    print("-" * 30)

    print("Step 2: Define the 'phase-shaping medium'.")
    # A phase-shaping medium applies a spatially-varying phase mask.
    # This is a SCALAR operation, meaning it multiplies both Ex and Ey
    # by the same complex number at each point (x,y).
    # Let's create a random phase mask T = exp(i * phi).
    random_phases = np.random.rand(2, 2) * 2 * np.pi
    T = np.exp(1j * random_phases)
    print("Transmission mask T (a scalar operator at each point):\n", T)
    print("-" * 30)

    print("Step 3: Apply the medium to the input beam.")
    # The output field E_out is T * E_in.
    # This operation is performed element-wise on each component.
    # E_out = [T * Ex_in, T * Ey_in]
    Ex_out = T * Ex_in
    Ey_out = T * Ey_in # Multiplying the mask by a zero field results in a zero field.
    
    # We output the numbers in the resulting components:
    print("Output Ex component (T * Ex_in):\n", Ex_out)
    print("Output Ey component (T * Ey_in):\n", Ey_out)
    print("-" * 30)

    print("Analysis:")
    print("The input beam had a zero Ey component.")
    # We check if the output Ey component is still zero.
    if np.all(Ey_out == 0):
        print("After passing through the scalar phase mask, the Ey component is still an array of zeros.")
        print("Free-space propagation is also a scalar operator and will likewise not create or mix polarization components.")
        print("\nConclusion: The described system is fundamentally unable to generate a y-polarized component from a purely x-polarized input.")
        print("Therefore, it cannot convert a uniformly polarized beam into a vector beam.")
    else:
        # This case is logically impossible given the definitions.
        print("An error occurred; the Ey component is unexpectedly non-zero.")

# Run the demonstration to explain the answer.
demonstrate_scalar_operation()