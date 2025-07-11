import math

def get_digits(x, num_digits):
    """
    Returns the digit sequence A for a number x.
    Note: This returns the standard representation, not dual ones like .999...
    """
    if not (0 <= x <= 10):
        raise ValueError("Input must be between 0 and 10")
    A = []
    # Special case for x=10
    if x == 10:
        return [9] * num_digits # Represents 9.999...
    
    a0 = int(x)
    A.append(a0)
    
    rem = x - a0
    for _ in range(num_digits - 1):
        rem *= 10
        digit = int(rem)
        A.append(digit)
        rem -= digit
    return A

def f_identity(x):
    """The identity function f(x) = x.
    This corresponds to T(A)_n = A_n / 10^n.
    Sum of T(A)_n is simply x."""
    return x

def f_digit_combo(x):
    """The function f(x) = 9*A_0 + A_1.
    This corresponds to T(A)_0 = 9*A_0 + A_1, and T(A)_n=0 for n > 0.
    The value of T(A)_n depends only on A_n, A_{n+1} as n=0.
    So it is shortsighted."""
    
    # We demonstrate this with two representations for an integer-adjacent number
    # like 3 = 2.999...
    if x == math.ceil(x) and x != 0:
        d = int(x - 1)
        # Sequence A' = (d+1, 0, 0, ...) -> (3, 0, 0, ...)
        val_A_prime = 9 * (d + 1) + 0
        
        # Sequence A = (d, 9, 9, ...) -> (2, 9, 9, ...)
        val_A = 9 * d + 9
        
        # We must check that they are equal for the function to be well-defined
        if val_A != val_A_prime:
            # This should not be reached for this function
            return float('nan')
        return val_A
    else:
        digits = get_digits(x, 3) # A_0, A_1 are sufficient
        A0 = digits[0]
        A1 = digits[1]
        return 9 * A0 + A1

def main():
    test_numbers = [math.pi, 0.5, 3, 9.9]
    print("For a given number x, we can calculate the values of different digitary functions.")
    print("f1(x) = x and f2(x) = 9*A_0 + A_1 are two linearly independent digitary functions.")
    print("We can construct infinitely many such functions, implying the dimension is countably infinite (N).")
    print("\n--- Example Calculations ---")
    for x in test_numbers:
        val1 = f_identity(x)
        val2 = f_digit_combo(x)
        digits = get_digits(x, 4)
        print(f"For x = {x} (digits begin {digits[0]}.{digits[1]}{digits[2]}{digits[3]}...):")
        print(f"  f1(x) = {val1}")
        print(f"  f2(x) = {val2} = 9*{digits[0]} + {digits[1]}")
    
    # Show an integer case explicitly for f2
    x = 4
    val2 = f_digit_combo(x)
    print(f"\nFor integer x = {x}:")
    print(f"  Representation {x}.00... gives 9*{x} + 0 = {9*x}")
    print(f"  Representation {x-1}.99... gives 9*{x-1} + 9 = {9*(x-1)+9}")
    print(f"  The values match, so f2 is a valid digitary function.")
    print(f"  f2({x}) = {val2}")

if __name__ == "__main__":
    main()