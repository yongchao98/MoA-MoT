import numpy as np

def decompose_field_vector():
    """
    This function demonstrates how to decompose a magnetic field fluctuation vector
    into components parallel and perpendicular to a local background magnetic field.

    This is the physically correct first step for calculating the helicity of a
    field-aligned wave like an AIC wave.
    """
    # Let's assume we are at L1 and have measured the average local magnetic field
    # in a standard coordinate system like GSE (where X is radial from Sun).
    # Note: This vector is NOT radial. It has strong Y and Z components,
    # representing a Parker Spiral orientation.
    # B_avg units are in nanoTeslas (nT).
    B_avg = np.array([-4.5, 3.8, -1.2])

    # Now, let's assume we measure a magnetic field fluctuation (the wave)
    # at a specific instant in the same coordinate system.
    delta_B = np.array([0.5, 0.8, 1.1])

    print(f"Local Average Magnetic Field B_avg (nT): {B_avg}")
    print(f"Measured Wave Fluctuation delta_B (nT): {delta_B}\n")

    # --- The Rigorous Method ---
    # 1. Find the direction of the local average magnetic field.
    #    This is the parallel direction.
    b_parallel_hat = B_avg / np.linalg.norm(B_avg)
    print(f"Direction of local field (unit vector): {np.round(b_parallel_hat, 3)}\n")

    # 2. Project the fluctuation vector (delta_B) onto the parallel direction
    #    to find the parallel component of the wave.
    #    This is done using the dot product.
    delta_B_parallel_scalar = np.dot(delta_B, b_parallel_hat)
    delta_B_parallel_vector = delta_B_parallel_scalar * b_parallel_hat
    
    # 3. The perpendicular component is what remains after subtracting the
    #    parallel component from the total fluctuation.
    #    This is the component used to calculate helicity.
    delta_B_perp_vector = delta_B - delta_B_parallel_vector
    
    print("--- Decomposition Results ---")
    print("For AIC wave analysis, we need components relative to the local field:")
    print(f"Parallel component of the wave (nT):   {np.round(delta_B_parallel_vector, 3)}")
    print(f"Perpendicular component of the wave (nT): {np.round(delta_B_perp_vector, 3)}")
    
    # The normalized magnetic helicity would then be calculated using a time series
    # of the components that make up delta_B_perp_vector, typically after
    # a Fourier transform. This vector demonstrates what components to use.
    # For example, one could define the two perpendicular axes and get the
    # components along them.

if __name__ == '__main__':
    decompose_field_vector()
