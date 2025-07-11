import textwrap

def explain_fuzzy_dimensions():
    """
    Illustrates the dimensional structure difference between Type-2 and Type-3
    fuzzy membership functions.
    """
    # --- Type-1 MF ---
    print("--- 1. Type-1 Membership Function (A 2D structure) ---")
    x1 = 5
    mu1 = 0.8
    print(textwrap.fill(
        "A Type-1 MF maps a 1D variable 'x' to a single 1D membership value 'mu'. This creates a simple 2D relationship, like a point on a graph.",
        width=80))
    print(f"Example: For x = {x1}, the membership is exactly mu = {mu1}")
    print("Dimensionality: (x, mu) -> 2D plot.\n")

    # --- Type-2 MF ---
    print("--- 2. Type-2 Membership Function (A 3D structure with a 2D uncertainty model) ---")
    x2 = 5
    mu2_lower = 0.6
    mu2_upper = 0.9
    print(textwrap.fill(
        "A Type-2 MF handles uncertainty about the membership value. For a single 'x', the membership is a range (an interval). This range is called the 'Footprint of Uncertainty' (FOU).",
        width=80))
    print(textwrap.fill(
        "While the full T2 MF is a 3D object, the MODEL of its uncertainty (the FOU) is a 2D region on the (x, mu) plane.",
        width=80))
    print(f"Example: For x = {x2}, the membership lies somewhere in the interval [{mu2_lower}, {mu2_upper}]")
    print("Dimensionality: The function is a 3D object. The uncertainty is modeled in 2D.\n")

    # --- Type-3 MF ---
    print("--- 3. Type-3 Membership Function (Adds 3D uncertainty modeling) ---")
    x3 = 5
    lower_bound_interval = (0.5, 0.7)
    upper_bound_interval = (0.8, 1.0)
    print(textwrap.fill(
        "A Type-3 MF adds another layer of uncertainty: the boundaries of the Type-2 FOU are now fuzzy themselves. The uncertainty model is no longer a flat 2D region but a volumetric, 3D region.",
        width=80))
    print(f"Example: For x = {x3}, the lower bound is the interval [{lower_bound_interval[0]}, {lower_bound_interval[1]}] and the upper bound is the interval [{upper_bound_interval[0]}, {upper_bound_interval[1]}]")
    print(textwrap.fill(
        "This gives 'thickness' or 'volume' to the FOU, effectively making the uncertainty representation itself three-dimensional.",
        width=80))
    print("Dimensionality Change: The uncertainty model transitions from 2D (a plane) to 3D (a volume).\n")
    
    # --- Conclusion ---
    print("--- CONCLUSION ---")
    print(textwrap.fill(
        "The fundamental difference is not just adding another variable, but that the structure used to model the uncertainty itself expands from a 2D plane to a 3D volume. Therefore, Type-3 introduces 'three-dimensional uncertainty modeling'.",
        width=80))

# Run the explanation
explain_fuzzy_dimensions()