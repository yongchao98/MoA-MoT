def find_ideal_ni_ce_ratio():
    """
    Summarizes findings on the ideal Ni/Ce ratio for Ni-Ceria catalysts
    in WGS and water splitting reactions based on scientific literature.
    """
    
    # The optimal range for Ni content is often cited in atomic percent (at%).
    # This is the percentage of Ni atoms relative to the total number of metal atoms (Ni + Ce).
    lower_bound_at_percent = 10
    upper_bound_at_percent = 20
    
    # A value around 15 at% is frequently cited as a "sweet spot" in many studies.
    representative_ni_at_percent = 15
    
    # Corresponding cerium atomic percent for the representative value.
    representative_ce_at_percent = 100 - representative_ni_at_percent
    
    # Calculate the atomic ratio (Ni/Ce) from the representative atomic percentages.
    # This represents the final equation requested.
    ni_ce_atomic_ratio = representative_ni_at_percent / representative_ce_at_percent
    
    print("Ideal Ni/Ce Ratio for Ni-Ceria Catalysts")
    print("-" * 45)
    print("The ideal ratio is not a single fixed value but falls within a well-studied range.")
    print(f"An optimal Ni content generally lies between {lower_bound_at_percent} at% and {upper_bound_at_percent} at% (atomic percent).")
    print("\nThis range maximizes the number of active Ni-Ceria interface sites while preventing Ni particle sintering and undesirable side reactions.")
    
    print("\n--- Example Calculation for a Representative Value ---")
    print(f"Let's consider a representative optimal value of {representative_ni_at_percent} at% Ni.")
    print(f"This corresponds to {representative_ce_at_percent} at% Ce.")
    print("\nThe Ni/Ce atomic ratio is calculated as:")
    print(f"Ratio = (Ni at%) / (Ce at%)")
    print(f"Ratio = {representative_ni_at_percent} / {representative_ce_at_percent}")
    print(f"Resulting Ni/Ce Atomic Ratio: {ni_ce_atomic_ratio:.3f}")
    
find_ideal_ni_ce_ratio()