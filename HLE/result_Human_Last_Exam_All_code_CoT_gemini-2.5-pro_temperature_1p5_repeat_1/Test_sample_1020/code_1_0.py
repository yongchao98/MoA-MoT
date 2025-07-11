import math

def solve():
    """
    Calculates the smallest complexity of a matrix product computing f_n,d.
    n, d are integers with 2 <= d <= n and d is odd.
    """
    try:
        # Taking input from the user
        n_str = input("Enter the integer n: ")
        d_str = input("Enter the odd integer d: ")
        n = int(n_str)
        d = int(d_str)

        # Validating the input
        if not (2 <= d <= n):
            print("Error: Input must satisfy 2 <= d <= n.")
            return
        if d % 2 == 0:
            print("Error: d must be an odd integer.")
            return

    except ValueError:
        print("Invalid input. Please enter integers only.")
        return

    # The formula for the minimal complexity C is:
    # C = 2 + 2n + (d-3)/2 * (floor(sqrt(n)) + ceil(n / floor(sqrt(n))))

    a = math.floor(math.sqrt(n))
    # Since n>=2, a >= 1, so no division by zero.
    b = math.ceil(n / a)

    # d is odd, so d-3 is even. Integer division // is appropriate.
    middle_sum = (d - 3) // 2 * (a + b)
    
    total_sum_m = 2 * n + middle_sum
    
    complexity = 2 + total_sum_m

    # Outputting the plan and the final equation step-by-step
    print("\nTo find the minimum complexity, we follow these steps:")
    print("1. Establish lower bounds on the intermediate matrix dimensions `m_i`.")
    print("   This gives m_1 >= n, m_{d-1} >= n, and m_i * m_{i+1} >= n.")
    print("2. Minimize the sum `m_1 + ... + m_{d-1}` subject to these constraints.")
    print("3. The optimal values are m_1 = n, m_{d-1} = n, and the intermediate dimensions m_2, ..., m_{d-2} alternate between a and b.")
    print(f"   a = floor(sqrt(n)) = floor(sqrt({n})) = {a}")
    print(f"   b = ceil(n / a) = ceil({n} / {a}) = {b}")
    
    print("\nThe total complexity is C = 2 + m_1 + ... + m_{d-1}")
    print(f"C = 2 + m_1 + m_{d-1} + (m_2 + ... + m_{d-2})")
    print(f"C = 2 + {n} + {n} + ((d-3)/2) * (a + b)")
    print(f"C = 2 + {2*n} + (({d}-3)/2) * ({a} + {b})")
    print(f"C = 2 + {2*n} + { (d-3)//2 } * {a+b}")
    print(f"C = 2 + {2*n} + {middle_sum}")
    print(f"C = {2 + 2*n + middle_sum}")
    print(f"C = {complexity}")
    
    print("\n---")
    print("The final expression for the minimum complexity is:")
    print(f"2 + 2*n + ((d-3)/2) * (floor(sqrt(n)) + ceil(n/floor(sqrt(n))))")
    print(f"\nFor n={n} and d={d}, the smallest complexity is:")
    print(f"<<<{complexity}>>>")

solve()
