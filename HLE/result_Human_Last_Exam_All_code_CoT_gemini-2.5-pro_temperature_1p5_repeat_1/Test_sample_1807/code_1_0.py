import sys

def calculate_fst_for_y_linked_marker():
    """
    Calculates the Fst value between males and females for a Y-linked marker.

    This demonstrates how sex determination systems can cause high genetic
    differentiation between sexes.

    Scenario:
    - We are analyzing a marker on the Y chromosome.
    - Males are XY, Females are XX.
    - The marker is present only in males.
    - We will assume the marker is fixed in the male population for simplicity.
    """
    print("Step 1: Define allele frequencies for a Y-linked marker in each sex.")
    # For a marker on the Y chromosome, its frequency in males is p_m.
    # Let's assume the marker is present in all males.
    p_m = 1.0
    print(f"The frequency of the marker in males (p_m) is: {p_m}")

    # Females (XX) do not have a Y chromosome, so the marker is absent.
    p_f = 0.0
    print(f"The frequency of the marker in females (p_f) is: {p_f}")
    print("-" * 30)

    print("Step 2: Calculate the average allele frequency (p_bar) across both sexes.")
    # p_bar is the mean of the subpopulation frequencies.
    p_bar = (p_m + p_f) / 2
    print(f"p_bar = (p_m + p_f) / 2")
    print(f"p_bar = ({p_m} + {p_f}) / 2 = {p_bar}")
    print("-" * 30)


    print("Step 3: Calculate the variance in allele frequency (Var(p)) between the sexes.")
    # Var(p) is the average of the squared differences from the mean frequency.
    var_p = ((p_m - p_bar)**2 + (p_f - p_bar)**2) / 2
    print(f"Var(p) = ((p_m - p_bar)^2 + (p_f - p_bar)^2) / 2")
    print(f"Var(p) = (({p_m} - {p_bar})^2 + ({p_f} - {p_bar})^2) / 2 = {var_p}")
    print("-" * 30)

    print("Step 4: Calculate Fst.")
    # Fst measures the proportion of total genetic variance contained in a subpopulation.
    # It is the ratio of the variance between subpopulations to the total variance.
    denominator = p_bar * (1 - p_bar)
    
    # Avoid division by zero if the marker is fixed or absent in the entire population.
    if denominator == 0:
      # If p_bar is 0 or 1, it means the allele is fixed or absent everywhere.
      # Fst is typically considered 0 in this case as there's no variation to partition.
      # However, in our case p_bar=0.5, so this check won't be triggered.
        fst = 0.0
    else:
        fst = var_p / denominator
        
    print(f"Fst = Var(p) / (p_bar * (1 - p_bar))")
    print(f"Fst = {var_p} / ({p_bar} * (1 - {p_bar}))")
    print(f"Fst = {var_p} / {denominator}")
    print("-" * 30)
    print(f"\nFinal calculated Fst value: {fst}")
    print("\nAn Fst of 1.0 indicates complete genetic differentiation between males and females for this marker.")
    print("This is a direct result of the marker being on a sex-specific chromosome (Y). A similar result would be found for W-linked markers in a ZW system. Therefore, the sex determination system is a key explanation.")


# Run the calculation
calculate_fst_for_y_linked_marker()
<<<B>>>