def calculate_crawley_nordstrom():
    """
    Calculates the Crawley-Nordstrom invariant for a given Calabi-Yau Link.
    """
    # The weights of the ambient space
    weights = [22, 29, 49, 50, 75]
    
    # The degree of the quasi-homogeneous polynomial
    # Based on the analysis, d = sum(weights) = 225
    d = 225
    
    # Calculate d^2
    d_sq = d**2
    
    # Calculate the sum of weights
    sum_w = sum(weights)
    
    # Calculate the sum of the squares of the weights
    sum_w_sq = sum(w**2 for w in weights)
    
    # Calculate the numerator and denominator of the invariant formula
    numerator = d_sq - sum_w_sq
    denominator = 3 * d - sum_w
    
    # The final invariant value
    invariant = numerator / denominator

    # --- Outputting the results as an equation ---
    
    # Create strings for the equation parts
    weights_sq_str = " + ".join([f"{w}^2" for w in weights])
    weights_sq_val_str = " + ".join([str(w**2) for w in weights])
    weights_sum_str = " + ".join(map(str, weights))

    print(f"Crawley-Nordström Invariant Formula: cn(X) = (d^2 - \u03A3w_i^2) / (3d - \u03A3w_i)")
    print("-" * 30)
    print(f"Given values:")
    print(f"d = {d}")
    print(f"w = {weights}")
    print("-" * 30)
    
    print("Step 1: Calculate the terms for the numerator")
    print(f"d^2 = {d}^2 = {d_sq}")
    print(f"\u03A3w_i^2 = {weights_sq_str}")
    print(f"      = {weights_sq_val_str}")
    print(f"      = {sum_w_sq}")
    print(f"Numerator = {d_sq} - {sum_w_sq} = {numerator}")
    print("-" * 30)

    print("Step 2: Calculate the terms for the denominator")
    print(f"3d = 3 * {d} = {3 * d}")
    print(f"\u03A3w_i = {weights_sum_str} = {sum_w}")
    print(f"Denominator = {3 * d} - {sum_w} = {denominator}")
    print("-" * 30)

    print("Step 3: Form the final equation and compute the result")
    # Using integer components for the final display
    print(f"cn(X) = ({d_sq} - {sum_w_sq}) / ({3 * d} - {sum_w})")
    print(f"      = {numerator} / {denominator}")
    print(f"      \u2248 {invariant}")
    print("-" * 30)
    print(f"The final Crawley-Nordström invariant is {numerator}/{denominator}")

calculate_crawley_nordstrom()

# The final answer in the required format
final_answer = 38774 / 450
print(f"\n<<<{final_answer}>>>")