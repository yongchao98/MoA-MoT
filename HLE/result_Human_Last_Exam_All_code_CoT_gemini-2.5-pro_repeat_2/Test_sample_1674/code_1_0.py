import numpy as np

def solve_optical_reversibility():
    """
    Analyzes the reversibility of an optical system with a birefringent layer.

    The theory of reversing an optical system to recover the initial input beam
    relies on the system's overall transmission operator being mathematically
    invertible. We can represent polarization-affecting components with 2x2
    Jones matrices. The reversibility of the system depends on whether these
    matrices have an inverse.
    """
    print("--- Analyzing Optical System Reversibility ---")

    # Case 1: The birefringent medium is a lossless component (e.g., a wave plate).
    # A wave plate shifts the phase between polarization components but doesn't
    # lose information. Its Jones matrix is invertible.
    print("\nCase 1: The added layer is a lossless wave plate.")
    
    # Jones matrix for a quarter-wave plate with its fast axis at 45 degrees.
    # The equation is: B = [[cos^2(a) + i*sin^2(a), cos(a)sin(a)(1-i)],
    #                      [cos(a)sin(a)(1-i), sin^2(a) + i*cos^2(a)]] for a=pi/4
    # which simplifies to:
    B_wave_plate = 0.5 * np.array([[1 - 1j, 1 + 1j],
                                  [1 + 1j, 1 - 1j]], dtype=complex)
    
    print("Jones Matrix for the Wave Plate (B_wave_plate):")
    for row in B_wave_plate:
        print(f"[{row[0]:.2f}, {row[1]:.2f}]")


    try:
        # Attempt to find the inverse matrix.
        B_wave_plate_inv = np.linalg.inv(B_wave_plate)
        print("\nResult: The wave plate matrix is INVERTIBLE.")
        print("The inverse represents the operation to reverse the wave plate's effect.")
        print("In this case, the theory of reversibility HOLDS.")
        
        print("\nInverse of B_wave_plate:")
        for row in B_wave_plate_inv:
            print(f"[{row[0]:.2f}, {row[1]:.2f}]")

    except np.linalg.LinAlgError:
        print("\nResult: The wave plate matrix is NOT invertible.")

    # Separator for clarity
    print("\n" + "="*60 + "\n")

    # Case 2: The birefringent medium is a lossy component (e.g., a polarizer).
    # A polarizer absorbs/removes one polarization component completely.
    # This is an irreversible loss of information. Its matrix is not invertible.
    print("Case 2: The added layer is a polarizer.")
    
    # Jones matrix for a perfect horizontal polarizer.
    # The equation is: B = [[1, 0], [0, 0]]
    B_polarizer = np.array([[1, 0],
                            [0, 0]], dtype=complex)

    print("Jones Matrix for the Polarizer (B_polarizer):")
    for row in B_polarizer:
        print(f"[{row[0]:.2f}, {row[1]:.2f}]")

    try:
        # Attempt to find the inverse matrix. This will fail.
        B_polarizer_inv = np.linalg.inv(B_polarizer)
        print("\nResult: The polarizer matrix is INVERTIBLE.")
    except np.linalg.LinAlgError:
        print("\nResult: The polarizer matrix is NOT INVERTIBLE (it is a singular matrix).")
        print("Because the polarizer deletes information (the vertical polarization component),")
        print("it's impossible to reverse the operation.")
        print("Therefore, when this layer is added, the theory of reversibility FAILS.")

solve_optical_reversibility()