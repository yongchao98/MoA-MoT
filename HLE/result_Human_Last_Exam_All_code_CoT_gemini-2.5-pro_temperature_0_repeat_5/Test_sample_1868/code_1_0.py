def calculate_fair_score():
    """
    Calculates the theoretical maximum FAIR compliance score R based on given parameters.
    """
    # Given parameters
    c = 0.95  # Consistency level of the identifier resolution
    b = 3     # Branching factor of semantic version control

    # --- Model the FAIR components ---

    # Findability (f) is limited by the consistency of the identifier resolution.
    # f_max = c
    f = c

    # Accessibility (a) is assumed to be ideal (1.0) due to "best implementation".
    a = 1.0

    # Interoperability (i) is inversely affected by the semantic branching factor.
    # A higher branching factor creates more ambiguity.
    # i_max = 1 / b
    i = 1 / b

    # Reusability (r) is assumed to be ideal (1.0) due to "best implementation".
    r = 1.0

    # --- Calculate the overall FAIR score R ---

    # R is the arithmetic mean of the four components.
    R = (f + a + i + r) / 4

    # --- Output the results ---
    print("The FAIR score R is calculated as the average of its components: R = (f + a + i + r) / 4")
    print(f"Based on the model:")
    print(f"  - Findability (f) = c = {f}")
    print(f"  - Accessibility (a) = Best Practice = {a}")
    print(f"  - Interoperability (i) = 1/b = 1/{b} = {i:.4f}")
    print(f"  - Reusability (r) = Best Practice = {r}")
    print("\nFinal Equation:")
    # The final print statement shows the equation with all numbers.
    print(f"R = ({f} + {a} + (1/{b}) + {r}) / 4")
    
    print(f"\nCalculated FAIR Score R = {R:.4f}")

calculate_fair_score()
<<<0.8208>>>