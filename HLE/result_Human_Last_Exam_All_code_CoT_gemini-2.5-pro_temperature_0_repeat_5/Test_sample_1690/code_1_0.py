import math

def solve_pioneer_problem():
    """
    Calculates the answers for the Pioneer probe problem based on the Wuxing architecture.
    """
    # Input parameters
    P = 100
    x = 4
    y = 5
    z = 1

    # a. Calculate camera shutdown time
    threshold_camera = x + y + z
    t_camera = 400 * (math.log(P) - math.log(threshold_camera)) / math.log(2)
    a = round(t_camera)

    # b. Calculate sensor shutdown time
    threshold_sensor = x + z
    t_sensor = 400 * (math.log(P) - math.log(threshold_sensor)) / math.log(2)
    b = round(t_sensor)

    # c. Calculate memory usage in D
    # 4 input ints (P,x,y,z) = 4 * 5D = 20D
    # 2 result ints (days_a, days_b) = 2 * 5D = 10D
    # 3 temporary fracs (for ln2, numerator, result) = 3 * 6D = 18D
    # Total = 20 + 10 + 18 = 48D
    c = 48

    # d. Calculate the number of exp function calls
    def count_exp_calls(v):
        """Calculates exp() calls for a given ln(v) based on the provided algorithm."""
        # The number of iterations is ceil(log2(v * 10^9))
        return math.ceil(math.log2(v * 10**9))

    # The program needs to compute ln for four values
    total_exp_calls = (count_exp_calls(2) +
                       count_exp_calls(P) +
                       count_exp_calls(threshold_camera) +
                       count_exp_calls(threshold_sensor))
    d = total_exp_calls
    
    # Print the final answer in the format a:b:c:d
    # This fulfills the requirement to "output each number in the final equation"
    print(f"{a}:{b}:{c}:{d}")

solve_pioneer_problem()