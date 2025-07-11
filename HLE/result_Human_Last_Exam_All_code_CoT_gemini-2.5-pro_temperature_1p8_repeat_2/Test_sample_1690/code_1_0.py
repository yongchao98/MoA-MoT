import math

def calculate_shutdown_times():
    # Given values
    P = 100.0
    x = 4.0
    y = 5.0
    z = 1.0

    # a. When to turn off the camera subsystem?
    # Power requirement for all systems
    C1 = x + y + z
    # The formula is t = 400 * log2(P/C)
    # Using python's math library: t1 = 400 * math.log2(P / C1)
    # With the Wuxing ln function: t1 = 400 * ln(P/C1) / ln(2)
    t1 = 400 * math.log(P / C1) / math.log(2)
    # Rounded to total days
    a = round(t1)

    # b. When to turn off the sensor subsystem?
    # Power requirement for control and sensor systems
    C2 = x + z
    # Using the same formula: t2 = 400 * log2(P/C2)
    t2 = 400 * math.log(P / C2) / math.log(2)
    # Rounded to total days
    b = round(t2)

    # c. What is the memory usage in D?
    # Data type sizes in Decimal positions (D)
    size_of_char = 2
    size_of_frac = 3 * size_of_char # (n, d, e are all char types)

    # Variable counts in a hypothetical C program
    # main() vars: P, x, y, z, C1, C2, num_400, val_2, val_P_C1, val_P_C2, ln_2, ln_P_C1, ln_P_C2, t1, t2 (15 frac)
    # Let's use a more streamlined count from a reasonable implementation:
    # P, x, y, z, C1, C2, const400, res_t1, res_t2, ln_2, ln_pc1, ln_pc2, v2, v_pc1, v_pc2 (15 fracs)
    # A lean implementation might have fewer. A more direct translation would use about 14 fracs.
    # e.g., P,x,y,z, C1,C2, val400, P_div_C1, P_div_C2, ln2, ln_PdivC1, ln_PdivC2, t1, t2 -> 14 fracs
    vars_main_frac = 14
    
    # ln(x) vars: l, r, m, y (4 frac)
    vars_ln_frac = 4

    # exp(x) vars: y, p (2 frac), n (1 char)
    vars_exp_frac = 2
    vars_exp_char = 1

    # Peak memory is when main calls ln, which calls exp
    mem_main = vars_main_frac * size_of_frac
    mem_ln = vars_ln_frac * size_of_frac
    mem_exp = vars_exp_frac * size_of_frac + vars_exp_char * size_of_char
    c = mem_main + mem_ln + mem_exp

    # d. What is the number of time this program call function exp?
    # exp() is called once per loop iteration in ln(x).
    # Number of iterations for ln(x) is ceil(log2(x / epsilon))
    epsilon = 1e-9
    
    # The program needs to compute ln(2), ln(P/C1), and ln(P/C2)
    # P/C1 = 100/10 = 10
    # P/C2 = 100/5 = 20
    
    calls_for_ln_2 = math.ceil(math.log2(2.0 / epsilon))
    calls_for_ln_10 = math.ceil(math.log2(10.0 / epsilon))
    calls_for_ln_20 = math.ceil(math.log2(20.0 / epsilon))

    d = calls_for_ln_2 + calls_for_ln_10 + calls_for_ln_20

    print(f"{int(a)}:{int(b)}:{int(c)}:{int(d)}")

calculate_shutdown_times()
#<<<1329:1729:122:100>>>