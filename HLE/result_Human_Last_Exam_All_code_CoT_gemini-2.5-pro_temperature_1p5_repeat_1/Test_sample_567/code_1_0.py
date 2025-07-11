def fibonacci(n):
    """
    Calculates the n-th Fibonacci number using an iterative method.
    Using the convention F_1=1, F_2=1, F_3=2, ...
    """
    if n <= 0:
        return 0
    elif n == 1:
        return 1
    a, b = 1, 1
    for _ in range(n - 2):
        a, b = b, a + b
    return b

def solve():
    """
    Calculates the value of 'a' where the volume constraint is the only
    obstruction for the symplectic embedding of the ellipsoid E(1,a) into a ball.
    """
    print("The condition that the only obstruction is volume means we are looking for 'a' such that c(a) = sqrt(a).")
    print("From symplectic geometry, this equality holds for a sequence a_k = (F_{2k+1} / F_{2k-1})^2, where F_n are Fibonacci numbers.")
    print("We will calculate the first non-trivial value, which corresponds to k=1.")
    print("-" * 20)

    k = 1
    
    # Calculate the indices for the Fibonacci numbers
    index_numerator = 2 * k + 1
    index_denominator = 2 * k - 1
    
    print(f"For k = {k}:")
    print(f"The index for the numerator is 2*k + 1 = {index_numerator}")
    print(f"The index for the denominator is 2*k - 1 = {index_denominator}")
    
    # Get the Fibonacci numbers
    f_num = fibonacci(index_numerator)
    f_den = fibonacci(index_denominator)
    
    print(f"The Fibonacci number F_{index_numerator} is {f_num}")
    print(f"The Fibonacci number F_{index_denominator} is {f_den}")

    # Calculate the ratio
    ratio = f_num / f_den
    
    # Calculate 'a'
    a = ratio ** 2
    
    print("\nThe final equation is:")
    # Using f-string to format the output as an equation
    print(f"a = (F_{index_numerator} / F_{index_denominator})^2 = ({f_num} / {f_den})^2 = {int(ratio)}^2 = {int(a)}")
    
    print("\nThe value of 'a' is:", int(a))

solve()