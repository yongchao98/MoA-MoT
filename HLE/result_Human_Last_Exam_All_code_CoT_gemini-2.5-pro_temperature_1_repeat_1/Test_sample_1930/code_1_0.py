import math

def get_digits(x, num_digits):
    """
    Computes the decimal digit sequence of a number x in [0, 10].
    A = (A_0, A_1, ..., A_{num_digits-1})
    x = A_0 + A_1/10 + A_2/100 + ...
    """
    if not (0 <= x <= 10):
        raise ValueError("Input x must be in [0, 10]")
    
    # Handle the edge case of x=10, represented as 9.999...
    if x == 10:
        return [9] * num_digits
        
    digits = []
    # A_0 is the integer part
    A_0 = math.floor(x)
    digits.append(A_0)
    
    remainder = x - A_0
    for _ in range(num_digits - 1):
        remainder *= 10
        A_n = math.floor(remainder)
        digits.append(A_n)
        remainder -= A_n
        
    return digits

def digitary_f_const(x, const=1, num_digits=100):
    """Implements f(x) = const using a shortsighted map."""
    A = get_digits(x, num_digits + 2) # Need A_n, A_{n+1}, A_{n+2}
    
    # g_0(d1, d2, d3) = const
    # g_n = 0 for n > 0
    
    # The sum is just the first term
    result = const
    
    print(f"For f(x) = {const}, x = {x:.4f}")
    print(f"  Digits A = {get_digits(x, 5)}...")
    print(f"  f(x) = g_0(A_0,A_1,A_2) = {const}")
    return result

def digitary_f_linear(x, c=1, num_digits=100):
    """Implements f(x) = c*x using a shortsighted map."""
    A = get_digits(x, num_digits + 2)
    
    # T(A)_n = g_n(A_n, A_{n+1}, A_{n+2}) = c * A_n / (10**n)
    
    total_sum = 0
    equation_str = []
    for n in range(num_digits):
        term = c * A[n] / (10**n)
        total_sum += term
        if n < 5: # Print first few terms of the equation
            equation_str.append(f"{c}*{A[n]}/10^{n}")

    print(f"For f(x) = {c}*x, x = {x:.4f}")
    print(f"  Digits A = {get_digits(x, 5)}...")
    print(f"  f(x) = {' + '.join(equation_str)} + ...")
    print(f"  Sum of first {num_digits} terms = {total_sum:.4f}")
    return total_sum
    
# Demonstrate for a few values
print("Demonstrating that f(x)=1 and f(x)=x are digitary functions.")
print("This implies the dimension is at least 2.")
print("The full analysis shows the dimension is countably infinite (N).")
print("-" * 20)
digitary_f_const(math.pi)
print("-" * 20)
digitary_f_linear(math.pi)
print("-" * 20)
# Test a value with dual representation
digitary_f_const(1.0)
print("-" * 20)
digitary_f_linear(1.0)
print("-" * 20)
# Test the other representation of 1.0 (0.999...)
# Our get_digits function doesn't produce this, so we'll simulate it.
print("f(x) = 1*x, for x = 0.999... = 1.0")
A_9s = [0] + [9]*100
total_sum_9s = 0
equation_str_9s = []
for n in range(100):
    term = 1 * A_9s[n] / (10**n)
    total_sum_9s += term
    if n < 5:
        equation_str_9s.append(f"1*{A_9s[n]}/10^{n}")
print(f"  Digits A = {A_9s[:5]}...")
print(f"  f(x) = {' + '.join(equation_str_9s)} + ...")
print(f"  Sum of first 100 terms = {total_sum_9s:.4f}")
