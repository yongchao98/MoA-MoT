import math

def calculate_cn(n):
    """
    Calculates the system-independent prefactor c_n for the fully
    f-connected Ree-Hoover diagram in the virial expansion.

    Args:
        n (int): The number of particles in the cluster (n >= 2).

    Returns:
        float: The value of the prefactor c_n.
    """
    if not isinstance(n, int) or n < 2:
        raise ValueError("n must be an integer of 2 or greater")
    
    if n == 2:
        # Special case for n=2
        return -1.0 / 2.0
    else:
        # General formula for n >= 3
        return ((-1.0)**n) / n

def main():
    """
    Main function to calculate and print c_n for a range of n.
    """
    print("The system-independent prefactor c_n for the fully f-connected Ree-Hoover diagram is:")
    
    # Calculate and print c_n for n from 2 to 10
    for n in range(2, 11):
        cn = calculate_cn(n)
        # The prompt asks to output each number in the final equation.
        # We interpret this as clearly presenting the result for each n.
        if n >= 3:
            numerator = f"(-1)^{n}"
            denominator = f"{n}"
            print(f"c_{n} = {numerator} / {denominator} = {(-1)**n} / {n} = {cn}")
        else: # n == 2
             print(f"c_{n} = -1 / 2 = {cn}")

if __name__ == "__main__":
    main()