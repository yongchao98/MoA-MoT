import math

def calculate_and_print_steps(M_str):
    """
    Calculates the optimal two-step gradient descent learning rates for a given M.

    Args:
        M_str (str): The condition number M as a string from user input.
    """
    try:
        kappa = float(M_str)
        if kappa < 1:
            print("Error: The condition number M must be a number greater than or equal to 1.")
            return
    except (ValueError, TypeError):
        print("Error: Please enter a valid number for M.")
        return

    # Print the general formulas
    print("\nThe best choice for the step sizes (gamma_1, gamma_2) is given by the formulas:")
    print("gamma_1 = (4*(M+1) - 2*sqrt(2)*(M-1)) / (M^2 + 6*M + 1)")
    print("gamma_2 = (4*(M+1) + 2*sqrt(2)*(M-1)) / (M^2 + 6*M + 1)")
    print("-" * 60)

    # Handle the trivial case M=1
    if kappa == 1.0:
        gamma1, gamma2 = 1.0, 1.0
        print("For M = 1, the problem is perfectly conditioned.")
        print("gamma_1 = 1.0")
        print("gamma_2 = 1.0")
        return
        
    print(f"Calculating for M = {kappa}:")
    
    # Breakdown of the calculation
    num_part1 = 4 * (kappa + 1)
    num_part2 = 2 * math.sqrt(2) * (kappa - 1)
    denominator = kappa**2 + 6*kappa + 1
    
    gamma1 = (num_part1 - num_part2) / denominator
    gamma2 = (num_part1 + num_part2) / denominator
    
    print("\n--- Calculation Breakdown ---")
    print(f"Numerator term (4*(M+1)):      4 * ({kappa}+1) = {num_part1}")
    print(f"Numerator term (2*sqrt(2)*(M-1)): 2 * {math.sqrt(2):.5f} * ({kappa}-1) = {num_part2:.5f}")
    print(f"Denominator (M^2 + 6*M + 1):    {kappa}^2 + 6*{kappa} + 1 = {denominator}")
    
    print("\n--- Final Equations with your value of M ---")
    print(f"gamma_1 = ({num_part1:.5f} - {num_part2:.5f}) / {denominator:.5f}")
    print(f"gamma_2 = ({num_part1:.5f} + {num_part2:.5f}) / {denominator:.5f}")

    print("\n--- Final Answer ---")
    print(f"gamma_1 = {gamma1}")
    print(f"gamma_2 = {gamma2}")


if __name__ == '__main__':
    # Since the problem is posed abstractly without a specific value for M,
    # the script prompts the user for one.
    try:
        m_input = input("Enter the condition number M (where M >= 1): ")
        calculate_and_print_steps(m_input)
    except (EOFError, KeyboardInterrupt):
        # Handle cases where input is not available (e.g., non-interactive script)
        # by running with a default example value.
        print("\nNo input provided. Running with an example value M=10.")
        calculate_and_print_steps('10')
