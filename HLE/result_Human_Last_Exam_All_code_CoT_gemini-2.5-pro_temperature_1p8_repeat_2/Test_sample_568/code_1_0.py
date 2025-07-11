import math

def calculate_mixing_repetitions(epsilon, delta):
    """
    Calculates the number of repetitions needed for two cups of liquid to be
    considered "the same" based on the given parameters.

    Args:
        epsilon (float): The maximum allowed difference in concentration (0 < epsilon < 1).
        delta (float): The fraction of liquid exchanged in each step (0 < delta < 1).
    """
    print(f"Given parameters: ε = {epsilon}, δ = {delta}")
    print("-" * 30)
    print("The condition for the liquids being 'the same' is |1 - 2δ|^n ≤ ε.")
    print("We solve for n, the number of repetitions.")

    # Handle the special case where delta is 0.5
    if delta == 0.5:
        n = 1
        print("\nSince δ = 0.5, the liquids mix perfectly in a single step.")
        print(f"Therefore, n = {n}")
        return n

    # The general case
    try:
        # Calculate the expression inside the ceiling function
        # n >= log(ε) / log(|1 - 2δ|)
        val = abs(1 - 2 * delta)
        log_eps = math.log(epsilon)
        log_val = math.log(val)
        result_float = log_eps / log_val
        n = math.ceil(result_float)
        
        print("\nThe formula for n is: ceil(log(ε) / log(|1 - 2δ|))")
        print("Substituting the values:")
        print(f"n = ceil(log({epsilon}) / log(|1 - 2 * {delta}|))")
        print(f"n = ceil({log_eps:.4f} / log({val:.4f}))")
        print(f"n = ceil({log_eps:.4f} / {log_val:.4f})")
        print(f"n = ceil({result_float:.4f})")
        print(f"n = {n}")

        return n

    except ValueError as e:
        print(f"\nError: Could not perform calculation. {e}")
        print("Please ensure epsilon and delta are within the (0, 1) range.")
        return None

if __name__ == '__main__':
    # Set the parameters for the problem
    # Epsilon (ε): must be between 0 and 1
    epsilon_param = 0.01
    # Delta (δ): must be between 0 and 1
    delta_param = 0.1

    calculate_mixing_repetitions(epsilon_param, delta_param)
