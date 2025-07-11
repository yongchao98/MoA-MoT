import math

def calculate_proton_position():
    """
    Calculates the position of a proton after three revolutions in
    perpendicular electric and magnetic fields.
    """
    # Given constants
    B_z = 5.0e-4  # Magnetic field in Tesla (T)
    E_x = 5.0     # Electric field in Volts/meter (V/m)
    q = 1.602e-19 # Proton charge in Coulombs (C)
    m = 1.673e-27 # Proton mass in kilograms (kg)
    # The initial velocity is not needed for calculating the position after an integer number of revolutions.

    # The motion is a drift in the y-direction. After 3 revolutions, x and z positions are 0.
    # The final y position is given by the formula: y = -(6 * pi * m * E_x) / (q * B_z^2)
    
    # Calculate the final y-position
    y_final = -(6 * math.pi * m * E_x) / (q * B_z**2)
    
    # The final position vector
    x_final = 0.0
    z_final = 0.0
    
    # Output the explanation and the final equation with values
    print("The final position of the proton is determined by the drift over three revolution periods.")
    print("The position vector is (x, y, z). Since the drift is along the y-axis, x and z are 0.")
    print("\nThe formula for the y-coordinate is:")
    print("y = - (6 * pi * m * E_x) / (q * B_z^2)")
    print("\nSubstituting the values into the equation:")
    # Using f-strings to format the output of the equation with numerical values
    print(f"y = - (6 * {math.pi} * {m} * {E_x}) / ({q} * ({B_z})^2)")
    
    # Calculate the numerator and denominator to show the calculation steps
    numerator = -(6 * math.pi * m * E_x)
    denominator = q * B_z**2
    print(f"y = {numerator} / {denominator}")
    
    print(f"\nFinal position (x, y, z) = ({x_final}, {y_final}, {z_final}) meters")

calculate_proton_position()