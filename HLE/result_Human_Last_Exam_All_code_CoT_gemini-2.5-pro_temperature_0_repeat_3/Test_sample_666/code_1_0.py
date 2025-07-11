import math

def count_poles_in_circle():
    """
    Counts the number of poles of the function f(z) inside the contour Gamma,
    which is assumed to be the circle |z|^2 = 200.

    The poles are at z = a + 2*pi*i*k, where a is an integer in [-2024, 2024]
    and k is any integer.

    The condition for a pole to be inside the circle is a^2 + (2*pi*k)^2 < 200.
    """
    R_squared = 200
    a_limit = 2024
    pole_count = 0

    # k = 0 case
    # a^2 < 200
    a_max_k0 = math.floor(math.sqrt(R_squared - 1e-9)) # Use a small epsilon for floating point issues
    # Number of integers a from -a_max to a_max
    # We must also respect the condition |a| <= 2024, which is satisfied since a_max < 2024
    count_k0 = 2 * a_max_k0 + 1
    pole_count += count_k0
    print(f"For k = 0, y = 0.00")
    print(f"  a^2 < {R_squared:.2f}")
    print(f"  |a| <= {a_max_k0}")
    print(f"  Number of poles: {count_k0}")

    # k > 0 cases (k < 0 are symmetric)
    k = 1
    while True:
        y_k = 2 * math.pi * k
        y_k_squared = y_k**2
        if y_k_squared >= R_squared:
            break
        
        a_squared_max = R_squared - y_k_squared
        a_max = math.floor(math.sqrt(a_squared_max))
        
        # Number of integers a for this k, respecting |a| <= 2024
        num_a = 2 * min(a_max, a_limit) + 1
        
        # Add counts for both k and -k
        pole_count += 2 * num_a
        
        print(f"For k = {k} and k = -{k}, y = +/-{y_k:.2f}")
        print(f"  a^2 < {a_squared_max:.2f}")
        print(f"  |a| <= {a_max}")
        print(f"  Number of poles: {2 * num_a}")

        k += 1
        
    print(f"\nTotal number of poles inside the contour: {pole_count}")
    
    # The integral is 2 * pi * i * pole_count.
    # The question likely asks for the number of poles.
    print(f"The value of the integral is 2 * pi * i * {pole_count} = {2*pole_count} * pi * i.")
    print(f"The question most likely asks for the number of poles, which is a real number.")
    print(f"Final Answer: {pole_count}")

count_poles_in_circle()