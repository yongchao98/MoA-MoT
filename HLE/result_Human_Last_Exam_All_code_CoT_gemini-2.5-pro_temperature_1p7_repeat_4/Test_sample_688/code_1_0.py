import math

def calculate_cn_prefactor(n):
    """
    Calculates the system-independent prefactor c_n for a given integer n.

    The prefactor c_n appears in the virial expansion contribution from the
    fully f-connected Ree-Hoover diagram. The formula is c_n = -(n-1)/n!.

    Args:
        n (int): The number of nodes (particles) in the diagram. 
                 The formula is physically meaningful for n >= 2.

    Returns:
        float: The numerical value of the prefactor c_n, or None if the
               input is invalid.
    """
    if not isinstance(n, int) or n < 2:
        print("Error: The virial coefficient B_n is defined for integers n >= 2.")
        return None
        
    numerator = n - 1
    denominator = math.factorial(n)
    
    # Output the explanation and calculation steps
    print(f"For n = {n}, the prefactor c_n is calculated as follows:")
    print(f"c_n = -(n - 1) / n!")
    print(f"c_{n} = -({n} - 1) / {n}!")
    print(f"c_{n} = -{numerator} / {denominator}")
    
    # Calculate the final result
    result = -numerator / denominator
    print(f"c_{n} = {result}")
    
    return result

if __name__ == '__main__':
    # You can change the value of n to calculate the prefactor for a different
    # number of particles. For example, n = 2, 3, 4, 5, etc.
    n_value = 5
    calculate_cn_prefactor(n_value)