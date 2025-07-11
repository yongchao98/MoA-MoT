import numpy as np

def check_polarization(E_field_H, E_field_V, tolerance=1e-9):
    """
    Checks if a beam is uniformly polarized or a vector beam.
    A simple check is to see if both orthogonal components have significant energy.
    """
    energy_H = np.sum(np.abs(E_field_H)**2)
    energy_V = np.sum(np.abs(E_field_V)**2)

    if energy_H > tolerance and energy_V > tolerance:
        return "Vector Beam"
    else:
        return "Uniformly Polarized Beam"

# 1. Define the parameters for our simulation grid.
grid_size = 128
x = np.linspace(-3, 3, grid_size)
X, Y = np.meshgrid(x, y)

# 2. Create the input beam.
# It is a Laguerre-Gaussian beam, but for simplicity of demonstrating the
# polarization effect, a simple Gaussian is sufficient.
# The key property is its uniform linear polarization.
# We represent it as a Jones vector [E_H, E_V].
amplitude_profile = np.exp(-(X**2 + Y**2))
E_H_in = amplitude_profile
E_V_in = np.zeros((grid_size, grid_size), dtype=complex)

print(f"Input beam type: {check_polarization(E_H_in, E_V_in)}")

# 3. Define the random phase-shaping medium (T) and its inverse (T_inv).
# A phase-shaping medium does not mix polarizations. We apply its effect
# to each component separately.
random_phase = 2 * np.pi * np.random.rand(grid_size, grid_size)
T = np.exp(1j * random_phase)
T_inv = np.conj(T) # The inverse of a pure phase transform is its conjugate.

# 4. Simulate the beam passing through the system.
# The problem describes: Input -> Prop -> T -> Prop -> T_inv -> Output.
# Free-space propagation does not mix polarizations. The key interaction is
# with the media T and T_inv. We apply their transformations component-wise.

# After passing through the first medium T:
E_H_intermediate = T * E_H_in
E_V_intermediate = T * E_V_in  # This remains zero.

# After passing through the second medium T_inv:
E_H_out = T_inv * E_H_intermediate
E_V_out = T_inv * E_V_intermediate # This also remains zero.

# Note that E_H_out simplifies to (T_inv * T * E_H_in), which is just E_H_in.
# The full system including propagation would result in a final beam that is
# a propagated version of the input, but its polarization state is unchanged.

print(f"Output beam type: {check_polarization(E_H_out, E_V_out)}")

print("\n--- Conclusion ---")
print("The analysis shows that a system composed of free-space propagation and phase-shaping media cannot convert a uniformly polarized beam into a vector beam.")
print("The reason is that none of these components can create a new polarization state from the existing one.")
print("The final answer to the question 'Can we get an arbitrary vector beam output?' is:")
print("No")
