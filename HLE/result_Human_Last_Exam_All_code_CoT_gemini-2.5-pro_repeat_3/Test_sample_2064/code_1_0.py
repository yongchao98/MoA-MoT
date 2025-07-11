import math

def l(a):
    """
    Calculates the value of l(a) for the case n=1.
    l(a) = - (2a / (1+a^2)) * (1 + 2*ln(1+a^2))
    """
    a_squared = a * a
    log_val = math.log(1 + a_squared)
    val = - (2 * a / (1 + a_squared)) * (1 + 2 * log_val)
    return val

def solve():
    """
    Calculates the floor of the sum of l(a_i) for the first 10 prime numbers.
    """
    # The first 10 prime numbers
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    
    total_sum = 0
    
    # We will build a string representation of the sum
    sum_str_parts = []
    
    print("Calculating l(a_i) for the first 10 prime numbers a_i:")
    for p in primes:
        val = l(p)
        total_sum += val
        print(f"l({p}) = {val}")
        sum_str_parts.append(f"({val:.3f})")

    print("\nThe full sum is:")
    sum_str = " + ".join(sum_str_parts)
    print(f"Sum = {sum_str} = {total_sum}")
    
    final_answer = math.floor(total_sum)
    print(f"\nThe floor of the sum is:")
    print(f"Floor({total_sum}) = {final_answer}")

solve()