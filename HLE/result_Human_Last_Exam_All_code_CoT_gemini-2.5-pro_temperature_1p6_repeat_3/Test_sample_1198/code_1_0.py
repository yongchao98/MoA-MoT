import decimal

# Use high precision for calculations
decimal.getcontext().prec = 50

def investigate_path(m, n_max, s1_choice):
    """
    Investigates if a real number 'a' can be found for floor(a^n) = n mod m
    by following a specific path of choices for the floor values s_n.

    Args:
        m (int): The modulus.
        n_max (int): How many steps to check.
        s1_choice (int): The starting choice for floor(a).
    """
    print(f"--- Investigating existence for modulo {m} ---")
    print(f"Trying path with initial choice: s_1 = floor(a) = {s1_choice}\n")

    if s1_choice % m != 1 % m:
        print(f"Error: Initial choice s_1={s1_choice} doesn't satisfy s_1 % {m} == 1 % {m}")
        return

    # Initialize the interval for 'a'
    # The interval is represented by its endpoints [low, high)
    low = decimal.Decimal(s1_choice)
    high = decimal.Decimal(s1_choice + 1)
    print(f"Step n=1: Let s_1 = {s1_choice}. This requires a in [{low}, {high}).")
    
    current_s = [s1_choice]

    for n in range(2, n_max + 1):
        # Calculate the required remainder for s_n = floor(a^n)
        required_rem = n % m
        
        # Determine the possible range for s_n
        s_n_low_bound = low ** n
        s_n_high_bound = high ** n
        
        # Find all possible integers in the range [s_n_low_bound, s_n_high_bound)
        # that satisfy the modulo condition.
        start_int = int(s_n_low_bound.to_integral_value(rounding=decimal.ROUND_CEILING))
        end_int = int(s_n_high_bound.to_integral_value(rounding=decimal.ROUND_FLOOR))
        # The range is exclusive on the right, so if high bound is an integer, exclude it.
        if s_n_high_bound == end_int:
             end_int -= 1

        possible_s_n = []
        for k in range(start_int, end_int + 1):
            if k % m == required_rem:
                possible_s_n.append(k)

        print(f"\nStep n={n}:")
        print(f"  a in [{float(low):.6f}, {float(high):.6f}) => a^{n} in [{float(s_n_low_bound):.6f}, {float(s_n_high_bound):.6f})")
        print(f"  Integers in this range for a^{n}: {list(range(start_int, end_int + 1))}")
        print(f"  We need s_{n} = floor(a^{n}) % {m} == {required_rem}")
        
        if not possible_s_n:
            print(f"\nFAILURE: No integer in the range [{start_int}, {end_int}] satisfies the condition.")
            print("This path of choices does not lead to a valid 'a'.")
            return
            
        # To continue this specific path, we must make a choice. Let's pick the smallest possible value.
        s_n_choice = possible_s_n[0]
        current_s.append(s_n_choice)
        print(f"  Possible s_{n} values: {possible_s_n}. We choose the smallest: s_{n} = {s_n_choice}")
        
        # Update the interval for 'a' based on the new constraint s_n_choice <= a^n < s_n_choice + 1
        n_inv = decimal.Decimal(1)/decimal.Decimal(n)
        new_low = decimal.Decimal(s_n_choice) ** n_inv
        new_high = (decimal.Decimal(s_n_choice) + 1) ** n_inv
        
        # Intersect with the previous interval
        low = max(low, new_low)
        high = min(high, new_high)

        if low >= high:
             print("\nFAILURE: The new interval for 'a' is empty.")
             return

        print(f"  The new interval for 'a' is [{float(low):.6f}, {float(high):.6f}).")

    print("\nSUCCESS: The path of choices was possible up to n_max.")

# Run the investigation for the modulo 3 case, with the smallest starting choice.
investigate_path(m=3, n_max=5, s1_choice=1)