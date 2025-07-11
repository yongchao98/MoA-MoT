def calculate_optimal_ni_ce_ratio():
    """
    Calculates and prints a commonly cited optimal Ni/Ce atomic ratio for catalysts.

    The ideal ratio for Ni-Ceria catalysts depends on many factors, but a
    composition of 10 atomic percent Ni is frequently reported in scientific
    literature as a sweet spot for maximizing performance in reactions like
    Water Gas Shift and water splitting.

    This function calculates the corresponding Ni-to-Ce atomic ratio.
    """

    # Define the composition in atomic percent (at%)
    ni_atomic_percent = 10
    ce_atomic_percent = 100 - ni_atomic_percent

    # Calculate the atomic ratio of Ni to Ce
    # Ratio = (atoms of Ni) / (atoms of Ce)
    ni_ce_atomic_ratio = ni_atomic_percent / ce_atomic_percent

    # --- Output ---
    print("Determining a commonly cited optimal Ni/Ce atomic ratio for catalytic performance.")
    print("-" * 70)
    print("A highly effective catalyst composition is often reported to contain 10 atomic percent Ni and 90 atomic percent Ce.")
    print("\nThe calculation for the Ni/Ce atomic ratio is as follows:")
    print("\nEquation: Ratio = (Atomic Percent of Ni) / (Atomic Percent of Ce)")
    print(f"Numbers:  Ratio = {ni_atomic_percent} / {ce_atomic_percent}")

    # Print the final calculated ratio, rounded to three decimal places for clarity
    print(f"\nThe resulting optimal Ni/Ce atomic ratio is approximately: {ni_ce_atomic_ratio:.3f}")
    print("-" * 70)
    print("This ratio promotes a high dispersion of active Ni sites on the Ceria support while minimizing particle sintering.")

if __name__ == '__main__':
    calculate_optimal_ni_ce_ratio()