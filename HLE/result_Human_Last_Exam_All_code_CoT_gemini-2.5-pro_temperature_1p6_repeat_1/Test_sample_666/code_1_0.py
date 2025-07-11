import math

def count_poles():
    """
    Counts the number of poles z = a + 2*pi*i*k inside the region x^2 + y^2 < 200.
    Poles are given by a in [-2024, 2024] and k integer.
    The condition is a^2 + (2*pi*k)^2 < 200.
    """
    num_poles = 0
    
    # Case k = 0
    max_a_sq_k0 = 200
    max_a_k0 = int(math.sqrt(max_a_sq_k0))
    # a ranges from -max_a_k0 to max_a_k0
    count_k0 = 2 * max_a_k0 + 1
    num_poles += count_k0
    
    print(f"For k = 0, condition is a^2 < 200.")
    print(f"sqrt(200) is approximately {math.sqrt(200):.2f}.")
    print(f"So, a is in [{-max_a_k0}, {max_a_k0}]. Number of poles: {count_k0}")
    
    # Case k != 0
    k = 1
    while True:
        k_term = (2 * math.pi * k)**2
        if k_term >= 200:
            break
        
        max_a_sq = 200 - k_term
        max_a = int(math.sqrt(max_a_sq))
        # Count for k and -k
        count_k = (2 * max_a + 1) * 2
        num_poles += count_k
        
        print(f"\nFor k = +/-{k}, condition is a^2 < 200 - (2*pi*{k})^2 = {max_a_sq:.2f}.")
        print(f"sqrt({max_a_sq:.2f}) is approximately {math.sqrt(max_a_sq):.2f}.")
        print(f"So, a is in [{-max_a}, {max_a}]. Number of poles for k=+/-%d: {count_k}" %k)
        
        k += 1
        
    return num_poles

total_poles = count_poles()
coefficient = total_poles * 2

print(f"\nTotal number of poles inside the contour is {total_poles}.")
print(f"The integral is 2 * pi * i * (number of poles).")
print(f"Value = 2 * pi * i * {total_poles} = {coefficient} * pi * i.")
print("\nThe final answer is:")
print(f"{coefficient} * pi * i")
