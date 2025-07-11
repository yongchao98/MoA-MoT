import math

def main():
    """
    This function solves the problem by first finding the dominant root of the
    characteristic equation of the recurrence, and then using it to compute the limit.
    """
    
    # The limit we need to evaluate simplifies to 2 * ln(r1), where r1 is the
    # dominant root of the characteristic equation r^4 - 88r - 57 = 0.
    # We use Newton's method to find this root.
    # Let f(r) = r^4 - 88r - 57. Then f'(r) = 4r^3 - 88.
    
    # We found a root exists between 4 and 5, so we can use 4.5 as an initial guess.
    r1 = 4.5
    for _ in range(10):  # 10 iterations provide sufficient precision
        r1 = r1 - (r1**4 - 88*r1 - 57) / (4*r1**3 - 88)

    # Now calculate the final value: 10^4 * 2 * ln(r1)
    val_10000 = 10**4
    val_2 = 2
    ln_r1 = math.log(r1)
    
    result = val_10000 * val_2 * ln_r1
    integer_part = math.floor(result)
    
    # Print the numbers involved in the calculation as requested
    print(f"The target expression is the integer part of: {val_10000} * {val_2} * ln(r1)")
    print(f"The dominant root found is r1 = {r1}")
    print(f"The expression evaluates to: {val_10000} * {val_2} * {ln_r1} = {result}")
    print(f"The integer part of the result is: {integer_part}")

if __name__ == "__main__":
    main()
