import math

def solve():
    """
    Calculates the answers for the Pioneer probe problem based on the Wuxing architecture.
    """
    # Given values
    P = 100.0
    x = 4.0
    y = 5.0
    z = 1.0
    half_life_days = 400.0

    # --- Part a: Camera subsystem shutdown time ---
    # Power required to run all three subsystems
    power_req_all = x + y + z
    # Using the formula: t = 400 * ln(P_initial / P_req) / ln(2)
    # This is equivalent to: t = 400 * log2(P_initial / P_req)
    t_a = half_life_days * math.log2(P / power_req_all)
    a = round(t_a)

    # --- Part b: Sensor subsystem shutdown time ---
    # Power required to run control and sensor subsystems
    power_req_sensor = x + z
    # Using the same formula
    t_b = half_life_days * math.log2(P / power_req_sensor)
    b = round(t_b)

    # --- Part c: Memory usage in D ---
    # Based on analysis of a memory-efficient C program on Wuxing.
    # Data type sizes: frac = 6D, char = 2D.
    # main() vars: 7 frac (P,x,y,z,ln2,temp,result) = 7 * 6D = 42D
    # ln() vars: 4 frac (l,r,m,y) = 4 * 6D = 24D
    # exp() vars: 2 frac (y,p) + 1 char (n) = 2*6D + 2D = 14D
    # Total max stack memory = 42 + 24 + 14
    c = 42 + 24 + 14

    # --- Part d: Number of exp() function calls ---
    # An efficient program calls ln() 3 times: ln(2), ln(10), and ln(20).
    # The number of exp() calls within ln(val) is ceil(log2(val * 1e9)).
    # log2(val * 1e9) = log2(val) + log2(1e9) = log2(val) + 9 * log2(10)
    log2_1e9 = 9 * math.log2(10)
    
    # Calls for ln(2)
    calls_ln2 = math.ceil(math.log2(2) + log2_1e9)
    # Calls for ln(10), where 10 = P / power_req_all
    calls_ln10 = math.ceil(math.log2(10) + log2_1e9)
    # Calls for ln(20), where 20 = P / power_req_sensor
    calls_ln20 = math.ceil(math.log2(20) + log2_1e9)
    
    d = calls_ln2 + calls_ln10 + calls_ln20

    # --- Final Output ---
    # The prompt requires printing each number in the final equation.
    # We will print the final calculated values for a, b, c, and d in the specified format.
    print(f"{a}:{b}:{c}:{d}")

solve()
<<<1329:1729:80:100>>>