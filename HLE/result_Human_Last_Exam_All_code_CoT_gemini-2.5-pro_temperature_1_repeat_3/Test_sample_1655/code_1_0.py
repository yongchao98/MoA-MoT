import numpy as np

def simulate_beam_restoration():
    """
    This simulation demonstrates that an optical system composed of a random medium (T)
    followed by its inverse (T_inv) restores the original input beam.
    """
    # 1. Define the spatial grid for our simulation
    grid_size = 128
    x = np.linspace(-1, 1, grid_size)
    y = np.linspace(-1, 1, grid_size)
    xx, yy = np.meshgrid(x, y)

    # 2. Create a simple input beam with a fixed polarization.
    # We only need to model the spatial amplitude and phase (the scalar part),
    # as the polarization is uniform and unaffected by phase-only elements.
    # Let's use a simple Gaussian beam as our input.
    w0 = 0.5  # Beam waist
    E_input = np.exp(-(xx**2 + yy**2) / w0**2)
    print("Step 1: An input beam with a fixed polarization is defined.")

    # 3. Create the random phase-shaping medium, T.
    # This medium applies a different phase shift at each point in space.
    # The transformation is represented by T = exp(i * phi(x,y)).
    random_phase = np.random.rand(grid_size, grid_size) * 2 * np.pi
    T = np.exp(1j * random_phase)
    print("Step 2: A random phase-shaping medium 'T' is created.")

    # 4. Propagate the beam through the first medium T.
    # The beam's phase is now completely scrambled.
    E_after_T = T * E_input
    print("Step 3: The input beam passes through medium 'T'.")

    # 5. Define the inverse medium, T_inv.
    # For a phase-only medium, the inverse is its complex conjugate.
    # T_inv * T = conj(T) * T = |T|^2 = 1 (Identity)
    T_inv = np.conj(T)
    print("Step 4: The inverse medium 'T_inv' is defined.")

    # 6. Propagate the scrambled beam through the inverse medium.
    E_output_2 = T_inv * E_after_T
    print("Step 5: The resulting beam passes through the inverse medium 'T_inv'.")

    # 7. Compare the final output with the original input.
    # We check if the system restored the original beam.
    # The "final equation" we are checking is: E_output_2 - E_input = 0
    # A good way to check this for arrays is to calculate the norm of the difference.
    # If the norm is zero (or very close to it), the arrays are identical.
    difference_norm = np.linalg.norm(E_output_2 - E_input)

    print("\n--- Simulation Result ---")
    print("The system is designed to make the final output identical to the input.")
    print("We verify this by calculating the norm of the difference between the two beams.")
    print("Final Equation: norm(E_output_2 - E_input) = Result")
    print(f"norm({E_output_2.__class__.__name__} - {E_input.__class__.__name__}) = {difference_norm:.15f}")

    if np.allclose(E_output_2, E_input):
        print("\nConclusion: The final output beam is a perfect replica of the input beam.")
        print("Since the input has a fixed linear polarization, the output must also have a fixed linear polarization.")
        print("A vector beam requires spatially varying polarization. Therefore, it is not possible to generate one with this system.")
    else:
        print("\nConclusion: The simulation failed to restore the beam.")

if __name__ == '__main__':
    simulate_beam_restoration()
