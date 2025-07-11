def calculate_minimum_gratings():
    """
    Calculates the minimum number of diffraction gratings for CTIS.

    The logic is based on the principles of computed tomography.
    """

    # In CTIS, to reconstruct the 2D spatial-spectral data for each point in the scene,
    # we need projections at various angles.
    # Projections along a single line (collinear) are insufficient.
    # We need projections along at least two independent axes to perform
    # a tomographic reconstruction.
    required_dispersion_axes = 2

    # A single 1D diffraction grating provides dispersion along one axis.
    gratings_per_axis = 1

    # The total minimum number of gratings is the product of the required
    # dispersion axes and the number of gratings needed per axis.
    minimum_gratings = required_dispersion_axes * gratings_per_axis

    print("The principle of computed tomography requires projections from multiple, non-collinear angles.")
    print("In CTIS, a 1D grating only provides collinear projections (along one axis).")
    print("To enable reconstruction, we must disperse light along a second, independent axis.")
    print("\nCalculation of the minimum required gratings:")
    print(f"Required dispersion axes for tomography = {required_dispersion_axes}")
    print(f"1D gratings needed per axis = {gratings_per_axis}")
    print(f"\nFinal Equation: {required_dispersion_axes} * {gratings_per_axis} = {minimum_gratings}")

calculate_minimum_gratings()