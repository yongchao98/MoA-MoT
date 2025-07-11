import math

def analyze_sugar_tube_experiment():
    """
    Analyzes the physics behind the D-glucose tube experiment and determines its appearance.
    """

    # --- 1. Define physical principles and illustrative parameters ---

    # The phenomenon is governed by Optical Rotatory Dispersion (ORD).
    # The angle of rotation 'alpha' is given by: alpha = [a] * L * c
    # where L is path length, c is concentration, and [a] is specific rotation.
    # [a] depends on wavelength (lambda). A simplified model is [a] is proportional to 1/lambda^2.

    # We can combine [a] and c into a single constant 'k' for our model.
    # Let's define an equation for the rotation angle in degrees:
    # rotation_angle = (k_constant / wavelength_nm**2) * path_length_m
    k_constant = 4000000  # An illustrative constant for effect (units: deg * nm^2 / m)
    path_length_m = 1.0     # The length of the tube

    # Define wavelengths for different colors.
    wavelength_red_nm = 650
    wavelength_blue_nm = 450

    # --- 2. Calculate the rotation for different colors ---

    # Equation for rotation of Red Light:
    # angle_red = (k_constant / wavelength_red_nm**2) * path_length_m
    angle_red = (k_constant / wavelength_red_nm**2) * path_length_m

    # Equation for rotation of Blue Light:
    # angle_blue = (k_constant / wavelength_blue_nm**2) * path_length_m
    angle_blue = (k_constant / wavelength_blue_nm**2) * path_length_m


    # --- 3. Print the analysis and conclusion ---

    print("Analysis of Light passing through D-Glucose Solution:\n")
    print("The appearance of the tube is determined by how the D-glucose solution interacts with different colors of light.")
    print("This is explained by Optical Rotatory Dispersion.\n")
    print("We can model the rotation angle with the equation: angle = (k / wavelength^2) * L")
    print("Using an illustrative constant k = {} and a tube length L = {} m:\n".format(k_constant, path_length_m))

    print("--- Rotation Calculation for Red Light ({} nm) ---".format(wavelength_red_nm))
    print("Final Rotation Angle = ({} / {}^2) * {}".format(k_constant, wavelength_red_nm, path_length_m))
    print("Result: Red light is rotated by {:.1f} degrees.".format(angle_red))
    print("-" * 50)

    print("--- Rotation Calculation for Blue Light ({} nm) ---".format(wavelength_blue_nm))
    print("Final Rotation Angle = ({} / {}^2) * {}".format(k_constant, wavelength_blue_nm, path_length_m))
    print("Result: Blue light is rotated by {:.1f} degrees.".format(angle_blue))
    print("-" * 50)

    print("\nConclusion:")
    print("Blue light's polarization rotates much more than red light's over the same distance.")
    print("At any point along the tube, an observer on the side sees the color whose polarization plane is optimally aligned for scattering.")
    print("Because this alignment angle changes with distance down the tube, the bands of color form helices.")
    print("This results in the visual appearance of a rainbow that progresses in a spiral.")


analyze_sugar_tube_experiment()