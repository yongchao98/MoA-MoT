import math

def check_mod3():
    """
    Checks the construction of 'a' for the mod 3 case and shows where it fails.
    """
    # n=1: floor(a) = 1 (mod 3). Let's pick floor(a) = 1.
    # a is in [1, 2)
    low = 1.0
    high = 2.0
    print(f"n=1: condition floor(a) = 1 (mod 3). Pick floor(a) = 1.")
    print(f"Interval for a: [{low}, {high})")

    # n=2: floor(a^2) = 2 (mod 3)
    # For a in [1, 2), a^2 is in [1, 4). floor(a^2) can be 1, 2, 3.
    # The only one == 2 (mod 3) is 2.
    # So floor(a^2) must be 2. a^2 must be in [2, 3).
    k = 2
    low = max(low, math.sqrt(k))
    high = min(high, math.sqrt(k + 1))
    print(f"\nn=2: condition floor(a^2) = 2 (mod 3). Must have floor(a^2) = {k}.")
    print(f"Interval for a: [{low}, {high})")

    # n=3: floor(a^3) = 0 (mod 3)
    # For a in [sqrt(2), sqrt(3)), a^3 is in [2.828, 5.196).
    # floor(a^3) can be 2, 3, 4, 5. Only 3 is 0 (mod 3).
    k = 3
    low = max(low, k**(1/3))
    high = min(high, (k + 1)**(1/3))
    print(f"\nn=3: condition floor(a^3) = 0 (mod 3). Must have floor(a^3) = {k}.")
    print(f"Interval for a: [{low}, {high})")

    # n=4: floor(a^4) = 1 (mod 3)
    # For a in [3^(1/3), 4^(1/3)), a^4 is in [3^(4/3), 4^(4/3)) = [4.32, 6.35).
    # floor(a^4) can be 4, 5, 6. Only 4 is 1 (mod 3).
    k = 4
    low = max(low, k**(1/4))
    high = min(high, (k + 1)**(1/4))
    print(f"\nn=4: condition floor(a^4) = 1 (mod 3). Must have floor(a^4) = {k}.")
    print(f"Interval for a: [{low}, {high})")

    # n=5: floor(a^5) = 2 (mod 3)
    # For a in the current interval, let's find the range for a^5.
    low_a5 = low**5
    high_a5 = high**5
    print(f"\nn=5: condition floor(a^5) = 2 (mod 3).")
    print(f"For a in [{low}, {high}), a^5 is in [{low_a5}, {high_a5}).")
    
    possible_k = [i for i in range(math.ceil(low_a5), math.ceil(high_a5)) if i % 3 == 2]
    
    print(f"Possible integer values for floor(a^5): from {math.ceil(low_a5)} to {math.floor(high_a5)}.")
    print(f"Integers in this range congruent to 2 (mod 3): {possible_k}")
    if not possible_k:
        print("No possible value for floor(a^5) satisfies the condition.")
        print("The construction fails.")

# To avoid running this in the final output, I will just output the conclusion.
# check_mod3()
