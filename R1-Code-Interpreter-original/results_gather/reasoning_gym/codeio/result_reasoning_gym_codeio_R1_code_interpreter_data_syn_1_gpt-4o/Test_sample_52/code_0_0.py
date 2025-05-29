import math

def calculate_digits(n):
    # Slow iterative method
    digits_slow = 0
    temp = abs(n)
    while True:
        temp = temp // 10
        digits_slow += 1
        if temp == 0:
            break
    
    # Logarithmic method
    digits_fast = 1 if n == 0 else math.floor(math.log10(abs(n)) + 1)
    
    # String conversion method
    digits_faster = len(str(abs(n)))
    
    return {
        "digits_slow": digits_slow,
        "digits_fast": digits_fast,
        "digits_faster": digits_faster
    }

# Test with a number having 20 digits
n = 10**19
result = calculate_digits(n)
print(result)