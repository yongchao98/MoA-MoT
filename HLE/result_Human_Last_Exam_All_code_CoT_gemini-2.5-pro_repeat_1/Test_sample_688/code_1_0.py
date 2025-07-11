import math

def double_factorial(n):
    """Computes the double factorial for odd integers n >= -1."""
    if n < -1 or n % 2 != 1:
        # According to the problem's context, we only need it for odd integers.
        # Let's define it as 0 for even numbers to avoid errors in intermediate steps
        # if the logic were different, though our formula 2n-5 is always odd.
        if n % 2 == 0: return 0 
        raise ValueError("Double factorial is typically defined for odd integers >= -1.")
    if n == -1:
        return 1
    # math.prod is available in Python 3.8+
    # For older versions, a loop would be needed.
    return math.prod(range(1, n + 1, 2))

def double_factorial_str(n):
    """Returns a string representation of the double factorial product."""
    if n < -1 or n % 2 != 1:
        return "N/A"
    if n == -1:
        return "1 [by definition]"
    if n == 1:
        return "1"
    return " * ".join(map(str, range(1, n + 1, 2)))

def calculate_and_print_c_n(n_values):
    """Calculates and prints the prefactor c_n for a list of n values."""
    print("The system-independent prefactor c_n is given by the formula:")
    print("c_n = - (2n - 5)!! / n\n")

    for n in n_values:
        if n < 2:
            print(f"n={n}: c_n is not defined for n < 2.\n")
            continue
            
        df_arg = 2 * n - 5
        df_val = double_factorial(df_arg)
        
        numerator = -df_val
        denominator = n
        
        print(f"--- Calculation for n = {n} ---")
        # Outputting each number in the equation as requested
        print(f"c_{n} = - (2 * {n} - 5)!! / {n}")
        print(f"   = - ({df_arg})!! / {n}")
        
        if df_arg > 1:
            print(f"   = - ({double_factorial_str(df_arg)}) / {n}")
        elif df_arg == -1:
             print(f"   = - ({double_factorial_str(df_arg)}) / {n}")

        print(f"   = {numerator}/{denominator}")
        
        # Simplify fraction for the final result
        common_divisor = math.gcd(abs(numerator), denominator)
        num_final = numerator // common_divisor
        den_final = denominator // common_divisor
        
        if den_final == 1:
            print(f"   = {num_final}\n")
        else:
            print(f"   = {num_final}/{den_final}\n")

# Execute the calculation for n from 2 to 7
calculate_and_print_c_n(range(2, 8))