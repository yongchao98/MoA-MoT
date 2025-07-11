import math

def calculate_critical_correlation():
    """
    This function calculates the critical amount of correlation required
    to balance potentiation and depression, leading to selectivity.

    The critical correlation rho_crit is given by the formula:
    rho_crit = 1 / sqrt(N_v * N_s)

    where N_v and N_s are the number of neurons in the input layers v and s.
    """
    # Define example parameters for the number of neurons in each layer.
    # These values can be changed to study different network sizes.
    N_v = 100
    N_s = 120

    # Calculate the critical correlation coefficient
    # The variance term sigma^2 cancels out in the final expression for the
    # correlation coefficient rho.
    rho_crit = 1 / math.sqrt(N_v * N_s)

    # Print the equation with the specific numbers used
    print(f"The critical correlation (rho_crit) is determined by the equation:")
    print(f"rho_crit = 1 / sqrt({N_v} * {N_s})")

    # Print the final calculated value
    print("\nCalculated Value:")
    print(f"rho_crit = {rho_crit}")

    # The final answer format as requested
    print(f"\n<<<final answer>>>")
    print(f"<<<{rho_crit}>>>")

if __name__ == "__main__":
    calculate_critical_correlation()