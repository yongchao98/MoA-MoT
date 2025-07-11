import math

def calculate_gromov_hausdorff_distance():
    """
    Calculates the Gromov-Hausdorff distance between the interval [0,1]
    and the unit circle with the intrinsic metric.
    """
    # Define the parameters of the two spaces
    # L: Length of the interval [0, 1]
    L = 1.0
    # C: Circumference of the unit circle (radius=1)
    C = 2 * math.pi

    print("Problem: Find the Gromov-Hausdorff distance between the interval [0, 1] and the unit circle.")
    print(f"  - The interval has length L = {L}")
    print(f"  - The unit circle has radius 1, so its circumference is C = 2 * pi ≈ {C:.4f}\n")

    # The formula for d_GH([0,L], S_C) for C >= 2L is:
    # d_GH = max(L/2, (C - 2L) / sqrt(8))
    print("Using the formula: d_GH = max(L/2, (C - 2*L) / sqrt(8))")
    print("-" * 50)
    
    # Calculate the two terms inside the max function
    term1 = L / 2
    term2 = (C - 2 * L) / math.sqrt(8)

    # Calculate the final distance
    gh_distance = max(term1, term2)
    
    # Print the breakdown of the final equation
    print("Calculating the first term:")
    print(f"  L / 2 = {L} / 2 = {term1}")
    
    print("\nCalculating the second term:")
    print(f"  (C - 2*L) / sqrt(8) = ({C:.4f} - 2*{L}) / {math.sqrt(8):.4f}")
    print(f"  = {C - 2*L:.4f} / {math.sqrt(8):.4f} = {term2:.4f}")

    print("\nFinding the maximum of the two terms:")
    print(f"  d_GH = max({term1:.4f}, {term2:.4f})")
    
    print("\nFinal Answer:")
    print(f"The Gromov-Hausdorff distance is {gh_distance:.4f}")

calculate_gromov_hausdorff_distance()