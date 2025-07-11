import math

def calculate_sum(n):
    """
    Calculates the closed form for the sum:
    S_n = sum_{k=0 to n} (2k+1)^5 * C(2k, k) * C(2(n-k), n-k)
    The closed form is: 4^n * (n+1) * (70n^3 - 65n^2 + 55n + 1)
    """
    if n < 0:
        return 0
    
    # Calculate the polynomial part
    poly_val = 70 * n**3 - 65 * n**2 + 55 * n + 1
    
    # Calculate the power of 4 part
    power_of_4 = 4**n
    
    result = power_of_4 * (n + 1) * poly_val
    
    print(f"For n = {n}:")
    print(f"The formula is 4^n * (n+1) * (70*n^3 - 65*n^2 + 55*n + 1)")
    print(f"S_{n} = 4^{n} * ({n}+1) * (70*({n})^3 - 65*({n})^2 + 55*({n}) + 1)")
    print(f"S_{n} = {power_of_4} * {n+1} * {poly_val}")
    print(f"The final result is: {result}")

def main():
    try:
        n_str = input("Enter an integer value for n (e.g., 2): ")
        n = int(n_str)
        calculate_sum(n)
    except ValueError:
        print("Invalid input. Please enter an integer.")

if __name__ == "__main__":
    main()
