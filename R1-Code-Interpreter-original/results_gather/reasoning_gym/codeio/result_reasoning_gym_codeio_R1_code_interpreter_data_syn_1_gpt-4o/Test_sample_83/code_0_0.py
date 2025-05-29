def numberOfFactors(n):
    counter = 0
    while n > 1:
        i = 2
        while i <= n:
            if n % i == 0:
                n /= i
                counter += 1
                break
            i += 1
    return counter

def liouville(n):
    if n == 1:
        return 1
    else:
        return (-1) ** numberOfFactors(n)

# Test the function with n = 1
print(liouville(1))