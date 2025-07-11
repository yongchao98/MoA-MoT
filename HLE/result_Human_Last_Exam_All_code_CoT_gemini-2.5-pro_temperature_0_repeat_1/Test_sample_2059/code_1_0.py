def solve():
    """
    Calculates the sum of the squares of the coefficients of the given polynomial expansion.
    The problem reduces to solving a linear recurrence relation C_n = 4*C_{n-1} + 2*c_{n-1}
    and c_n = 3*C_{n-1} + 2*c_{n-1}, where C_n is the sum of squares for a product of n terms.
    The solution can be expressed as C_n = 4*A_{n-1} + 10*B_{n-1}, where (3+sqrt(7))^(n-1) = A_{n-1} + B_{n-1}*sqrt(7).
    We need to find C_20, which requires computing A_19 and B_19.
    """

    # Let (3+sqrt(7))^n = A_n + B_n * sqrt(7).
    # We start with n=1.
    # A_1 = 3, B_1 = 1.
    A = 3
    B = 1

    # We need A_19 and B_19. We loop from n=2 to 19.
    # This is a total of 18 iterations.
    for _ in range(18):
        # Recurrence relations:
        # A_{n+1} = 3*A_n + 7*B_n
        # B_{n+1} = A_n + 3*B_n
        A_next = 3 * A + 7 * B
        B_next = A + 3 * B
        A = A_next
        B = B_next

    # After the loop, A is A_19 and B is B_19.
    A_19 = A
    B_19 = B

    # The sum of squares is C_20 = 4 * A_19 + 10 * B_19.
    result = 4 * A_19 + 10 * B_19

    print(f"The sum of squares is given by C_20 = 4 * A_19 + 10 * B_19.")
    print(f"where (3+sqrt(7))^19 = A_19 + B_19 * sqrt(7).")
    print(f"The calculated values are:")
    print(f"A_19 = {A_19}")
    print(f"B_19 = {B_19}")
    print(f"The final equation is: 4 * {A_19} + 10 * {B_19} = {result}")

solve()