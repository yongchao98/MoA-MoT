def calculate_minimum_gratings():
    """
    Calculates the minimum number of diffraction gratings for single-shot
    computed tomography imaging spectrometry (CTIS).

    In CTIS, a 3D data cube (x, y, lambda) is reconstructed from a 2D sensor image.
    This tomographic reconstruction requires projections with angular diversity.
    """

    # For a full tomographic reconstruction of a 2D spatial scene, we need
    # to create projections that are dispersed along two independent axes.
    required_dispersion_axes = 2

    # A single, standard diffraction grating provides dispersion along only one axis.
    axes_per_grating = 1

    # The minimum number of gratings is the total required axes divided by the
    # number of axes each grating provides.
    min_gratings = required_dispersion_axes // axes_per_grating

    print("--- CTIS Grating Requirement Calculation ---")
    print(f"To reconstruct a spectral volume from a single 2D image, projections must be dispersed across a 2D plane.")
    print(f"Number of independent dispersion axes required: {required_dispersion_axes}")
    print(f"Number of dispersion axes provided by a single standard grating: {axes_per_grating}")
    print(f"\nFinal Equation: {required_dispersion_axes} (required axes) / {axes_per_grating} (axes per grating) = {min_gratings} gratings")
    print(f"\nTherefore, the minimum number of diffraction gratings necessary is {min_gratings}.")

calculate_minimum_gratings()