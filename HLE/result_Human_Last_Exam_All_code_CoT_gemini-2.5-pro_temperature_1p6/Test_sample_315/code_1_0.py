import math

def calculate_optical_rotation():
    """
    Calculates and displays the optical rotation for different wavelengths of light
    in a D-glucose solution, explaining the rainbow spiral phenomenon.
    """
    # Parameters of the experiment
    # Path length in decimeters (1 m = 10 dm)
    path_length_dm = 10.0
    # An example concentration of D-glucose solution in g/mL
    # (e.g., 200g in 1L of water = 0.2 g/mL)
    concentration_g_ml = 0.2

    # Specific rotation [alpha] in degrees*mL/(g*dm) for D-glucose at different wavelengths.
    # This value is an intrinsic property of the substance.
    # Note: Blue light is rotated more than red light.
    specific_rotations = {
        "Red (650 nm)": 45.0,
        "Yellow (589 nm)": 52.7,
        "Blue (436 nm)": 102.8
    }

    print("The formula for optical rotation is: Angle = [alpha] * length * concentration\n")
    print(f"Assuming a D-glucose solution with:")
    print(f"  - Path Length: {path_length_dm} dm (which is 1m)")
    print(f"  - Concentration: {concentration_g_ml} g/mL\n")

    print("--- Calculated Total Rotation at the end of the tube ---")

    for color, alpha in specific_rotations.items():
        # Calculate the total rotation angle using the formula
        total_rotation_deg = alpha * path_length_dm * concentration_g_ml
        
        # Calculate the number of full 360-degree turns
        num_rotations = total_rotation_deg / 360.0

        print(f"\nFor {color}:")
        print(f"  Angle = {alpha} * {path_length_dm} * {concentration_g_ml} = {total_rotation_deg:.1f} degrees")
        print(f"  This corresponds to {num_rotations:.2f} full spiral turns along the tube.")
        
    print("\nThis calculation shows that different colors rotate by vastly different amounts.")
    print("This difference in rotation for each color is what creates the spiral rainbow effect.")

if __name__ == "__main__":
    calculate_optical_rotation()