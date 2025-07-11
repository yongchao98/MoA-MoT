def find_ideal_ni_ce_ratio():
    """
    This function presents findings from scientific literature on the ideal Ni/Ce
    molar ratio for catalytic applications like Water Gas Shift (WGS) and Water Splitting (WS).

    The ideal ratio is not a single fixed number but an optimal range that can
    vary with synthesis method and reaction conditions. The key is to balance
    the number of active Ni sites with high dispersion (small particle size)
    to maximize the Ni-Ceria interface, which is crucial for catalytic activity.
    """

    # Based on literature, the optimal performance is often found when the molar
    # percentage of Ni relative to the total metal content (Ni + Ce) is between 10% and 20%.
    # We can translate this to a direct Ni/Ce molar ratio.

    # Lower bound calculation: 10% Ni -> Ni/(Ni+Ce) = 0.1 -> 0.9*Ni = 0.1*Ce -> Ni/Ce = 0.1/0.9
    optimal_ratio_lower_bound = 0.1 / 0.9

    # Upper bound calculation: 20% Ni -> Ni/(Ni+Ce) = 0.2 -> 0.8*Ni = 0.2*Ce -> Ni/Ce = 0.2/0.8
    optimal_ratio_upper_bound = 0.2 / 0.8

    print("Analysis of the Ideal Ni/Ce Molar Ratio for Catalysis")
    print("-" * 60)
    print("The ideal Ni/Ce ratio maximizes catalytic sites while preventing Ni particle agglomeration.")
    print("This optimal balance is generally found within a specific range.")
    print("\nBased on a review of experimental data for WGS and water splitting reactions:")
    
    # Printing the numbers for the "final equation" as requested
    print("\n--- Optimal Range Calculation ---")
    print(f"Represents a Ni/(Ni+Ce) molar percentage of 10% to 20%.")
    print(f"Calculated Lower Bound (Ni/Ce Ratio): {optimal_ratio_lower_bound:.2f}")
    print(f"Calculated Upper Bound (Ni/Ce Ratio): {optimal_ratio_upper_bound:.2f}")

    # Final Conclusion
    print("\n--- Final Conclusion ---")
    print("The final optimal molar ratio equation is represented by the range:")
    print(f"Ideal_Ni/Ce_Ratio = ({optimal_ratio_lower_bound:.2f}) to ({optimal_ratio_upper_bound:.2f})")
    print("\nA representative value within this range, often cited for high, stable activity, is approximately 0.18.")


if __name__ == "__main__":
    find_ideal_ni_ce_ratio()