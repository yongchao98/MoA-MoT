import numpy as np
from scipy.spatial.transform import Rotation

def solve_rotation_convention():
    """
    Finds the Euler angle convention that matches a given Tait-Bryan rotation.
    """
    # 1. Define the initial Tait-Bryan rotation and compute its matrix.
    # Convention: Extrinsic 'XYZ', Angles: alpha=10, beta=10, gamma=10 degrees.
    tb_angles = [10, 10, 10]
    tb_convention = 'XYZ'
    
    # SciPy's uppercase conventions are extrinsic.
    r_ref = Rotation.from_euler(tb_convention, tb_angles, degrees=True)
    m_ref = r_ref.as_matrix()

    print("Step 1: Calculate the reference rotation matrix from the Tait-Bryan angles.")
    print(f"Convention: Extrinsic {tb_convention}, Angles: α={tb_angles[0]}°, β={tb_angles[1]}°, γ={tb_angles[2]}°")
    print("Reference Rotation Matrix (R_ref):")
    print(np.round(m_ref, 4))
    print("-" * 50)

    # 2. Define the candidate Euler angles and conventions.
    euler_angles = [139.13, 14.11, -141.05]
    
    conventions_to_test = {
        'A': 'XZX',
        'B': 'XYZ',
        'C': 'YXY',
        'D': 'YZY',
        'E': 'ZYZ',
        'F': 'ZXZ'
    }

    print("Step 2 & 3: Test each convention with the given Euler angles.")
    print(f"Euler Angles: α'={euler_angles[0]}°, β'={euler_angles[1]}°, γ'={euler_angles[2]}°\n")

    correct_choice = None
    # 3 & 4. Iterate, calculate test matrices, and compare.
    for choice, conv_str in conventions_to_test.items():
        r_test = Rotation.from_euler(conv_str, euler_angles, degrees=True)
        m_test = r_test.as_matrix()

        print(f"Testing Choice {choice}: Convention '{conv_str}'")
        print("Resulting Matrix (R_test):")
        print(np.round(m_test, 4))

        # Compare the test matrix with the reference matrix
        is_match = np.allclose(m_ref, m_test, atol=1e-4)
        print(f"Comparison: Is R_test numerically equal to R_ref? -> {is_match}")
        
        if is_match:
            correct_choice = choice
        print("-" * 50)

    if correct_choice:
        print(f"Conclusion: The correct convention is {correct_choice} ({conventions_to_test[correct_choice]}).")
        print(f"The rotation from ZYZ(α'={euler_angles[0]}°, β'={euler_angles[1]}°, γ'={euler_angles[2]}°) is equivalent to XYZ(α={tb_angles[0]}°, β={tb_angles[1]}°, γ={tb_angles[2]}°).")
    else:
        print("Conclusion: No matching convention was found among the choices.")

solve_rotation_convention()