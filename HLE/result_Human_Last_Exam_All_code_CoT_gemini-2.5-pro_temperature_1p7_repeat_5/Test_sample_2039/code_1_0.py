import numpy as np

def get_rotation_matrix(axis, angle_deg):
    """
    Returns the 3x3 rotation matrix for a given axis and angle.
    - axis: 'x', 'y', or 'z'
    - angle_deg: angle in degrees
    """
    angle_rad = np.deg2rad(angle_deg)
    c, s = np.cos(angle_rad), np.sin(angle_rad)
    
    if axis == 'x':
        return np.array([
            [1, 0, 0],
            [0, c, -s],
            [0, s, c]
        ])
    elif axis == 'y':
        return np.array([
            [c, 0, s],
            [0, 1, 0],
            [-s, 0, c]
        ])
    elif axis == 'z':
        return np.array([
            [c, -s, 0],
            [s, c, 0],
            [0, 0, 1]
        ])
    else:
        raise ValueError("Axis must be 'x', 'y', or 'z'")

def main():
    # Step 1: Calculate the source rotation matrix
    # Extrinsic X-Y-Z rotation with alpha=beta=gamma=10 degrees.
    # Matrix multiplication order is Rz @ Ry @ Rx.
    alpha, beta, gamma = 10.0, 10.0, 10.0
    Rx_source = get_rotation_matrix('x', alpha)
    Ry_source = get_rotation_matrix('y', beta)
    Rz_source = get_rotation_matrix('z', gamma)
    R_source = Rz_source @ Ry_source @ Rx_source

    print("--- Task: Find Equivalent Euler Angle Convention ---")
    print(f"Initial extrinsic rotation: X({alpha} deg), Y({beta} deg), Z({gamma} deg)")
    
    # Step 2: Define target angles and conventions to test
    ap, bp, gp = 139.13, 14.11, -141.05
    print(f"Target proper Euler angles: alpha'={ap}, beta'={bp}, gamma'={gp}\n")

    # The order of angles (ap, bp, gp) corresponds to the order of axes in the convention.
    # E.g., for ZXZ: first rotation is around Z by ap, second around X by bp, third around Z by gp.
    conventions_to_test = {
        'A': ('x', 'z', 'x'),
        'B': ('x', 'y', 'z'),
        'C': ('y', 'x', 'y'),
        'D': ('y', 'z', 'y'),
        'E': ('z', 'y', 'z'),
        'F': ('z', 'x', 'z')
    }
    
    found_convention = None
    
    # Step 3: Test each convention
    for key, axes in conventions_to_test.items():
        axis1, axis2, axis3 = axes
        
        # Calculate R_target assuming an intrinsic convention
        R1 = get_rotation_matrix(axis1, ap)
        R2 = get_rotation_matrix(axis2, bp)
        R3 = get_rotation_matrix(axis3, gp)
        R_target = R1 @ R2 @ R3
        
        # Check if the calculated matrix matches the source matrix
        if np.allclose(R_source, R_target, atol=1e-4):
            found_convention = f"{key}. {''.join(axes).upper()}"
            R_matched = R_target
            break
            
    if found_convention:
        print("Found a match!")
        print(f"The correct convention is: {found_convention}\n")
        print("Comparing the rotation matrices:\n")
        print("Source Matrix from X_10 Y_10 Z_10 (extrinsic):")
        print(np.round(R_source, 4))
        print(f"\nTarget Matrix from {found_convention.split('. ')[1]}({ap}, {bp}, {gp}) (intrinsic):")
        print(np.round(R_matched, 4))
    else:
        print("Could not find a matching convention among the options.")

if __name__ == "__main__":
    main()