def calculate_four_terminal_conductance():
    """
    This script calculates and presents the formula for the four-terminal
    conductance G_12,34 of the specified quantum Hall device.
    """
    # Define symbolic representations for the physical quantities
    M = "M"  # Total number of spin-degenerate edge states
    N = "N"  # Number of reflected edge states by the QPC
    G0 = "e^2/h" # The conductance quantum

    print("Calculation of the four-terminal conductance G_12,34")
    print("-" * 50)
    print("The calculation is based on the Landauer-Büttiker formalism for quantum Hall systems.")
    print("Key intermediate results:")
    print(f"1. The relationship between the source current I (from 1 to 2) and terminal voltages V1, V2:")
    print(f"   I = ({G0}) * ({M} - {N}) * (V1 - V2)")
    print()
    print(f"2. The measured voltage difference between terminals 3 and 4:")
    print(f"   V3 - V4 = ({N} / {M}) * (V1 - V2)")
    print("-" * 50)
    print("The conductance G_12,34 is the ratio I / (V3 - V4).")
    print("Substituting the expressions and simplifying gives the final equation:")

    # Construct and print the final formula, highlighting each component
    numerator = f"{M} * ({M} - {N})"
    denominator = f"{N}"
    
    print(f"\nG_12,34 = ({G0}) * ({numerator}) / ({denominator})")
    
    print("\nWhere the terms in the final equation represent:")
    print(f"  - {G0} : The fundamental conductance quantum.")
    print(f"  - {M}          : The total number of edge states.")
    print(f"  - {M} - {N}      : The number of transmitted edge states through the QPC.")
    print(f"  - {N}          : The number of reflected edge states at the QPC.")

# Run the calculation and display the results
calculate_four_terminal_conductance()