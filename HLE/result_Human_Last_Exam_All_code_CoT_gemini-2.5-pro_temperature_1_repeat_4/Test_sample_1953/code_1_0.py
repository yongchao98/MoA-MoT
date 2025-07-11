import cmath

def solve_quantum_sensing():
    """
    Calculates the difference between 1 and the Quantum Fisher Information (QFI)
    for a distributed quantum sensing scenario.
    """
    try:
        d_str = input("Enter the number of sensor nodes (d): ")
        d = int(d_str)
        if d <= 0:
            print("Error: The number of sensor nodes (d) must be a positive integer.")
            return

        f_str = input("Enter the fidelity of the noisy GHZ state (F): ")
        F = float(f_str)
        if not (0.0 <= F <= 1.0):
            print("Error: The fidelity (F) must be between 0 and 1.")
            return

    except ValueError:
        print("Error: Invalid input. Please enter a positive integer for d and a number between 0 and 1 for F.")
        return

    # The Quantum Fisher Information (QFI) for the parameter theta is H = 4 * d * (2*F - 1)^2.
    # Here we show the calculation step-by-step.
    
    print("\nStep 1: Calculate the term (2*F - 1)")
    term1_val = 2 * F - 1
    print(f"(2 * {F} - 1) = {term1_val}")

    print("\nStep 2: Square the result")
    term2_val = term1_val ** 2
    print(f"({term1_val})^2 = {term2_val}")

    print("\nStep 3: Calculate the Quantum Fisher Information (QFI), H = 4 * d * (2*F - 1)^2")
    qfi = 4 * d * term2_val
    print(f"H = 4 * {d} * {term2_val} = {qfi}")
    
    # The final result is the difference between 1 and the QFI.
    result = 1 - qfi

    print("\nFinal Step: Calculate the difference between 1 and the QFI")
    print(f"Result = 1 - H")
    print(f"Result = 1 - {qfi} = {result}")

# Execute the function
solve_quantum_sensing()