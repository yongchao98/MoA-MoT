def calculate_minimal_conformal_weight():
    """
    Calculates the minimal non-zero conformal weight for V(p) at p=2
    based on the provided decomposition structure.
    """
    p = 2
    # The minimal non-zero conformal weight corresponds to the n=1 sector.
    n = 1

    print("Calculation of the minimal non-zero conformal weight for p=2:")
    print("The formula for the conformal weight of the highest-weight state in the n-th sector is h_n = (p * n * (n+2)) / 4.")
    print(f"For the minimal non-zero weight, we use p = {p} and n = {n}.")

    # Perform the calculation
    numerator = p * n * (n + 2)
    denominator = 4
    minimal_weight = numerator / denominator

    # Output the steps of the calculation
    print(f"h_1 = ({p} * {n} * ({n} + 2)) / {denominator}")
    print(f"h_1 = ({p} * {n+2}) / {denominator}")
    print(f"h_1 = {numerator} / {denominator}")
    print(f"h_1 = {minimal_weight}")

calculate_minimal_conformal_weight()