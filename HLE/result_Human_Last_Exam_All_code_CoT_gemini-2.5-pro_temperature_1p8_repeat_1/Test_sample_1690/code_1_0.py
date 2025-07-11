import math

def solve_pioneer_problem():
    """
    Calculates the answers for the Pioneer probe problem based on the Wuxing architecture.
    """
    # Given values
    P = 100
    x = 4
    y = 5
    z = 1

    print("Step-by-step calculation:")

    # --- Part a: Camera shutdown time ---
    power_all = x + y + z
    ratio_cam = power_all / P
    # Using the formula: t = -400 * ln(ratio) / ln(2)
    # Since ratio_cam = 0.1, we calculate -ln(1/0.1) = -ln(10)
    t_cam = -400 * math.log(ratio_cam) / math.log(2)
    a = round(t_cam)
    print(f"\na. The camera is turned off when power drops below {power_all}W.")
    print(f"   The equation is t = -400 * ln(({x} + {y} + {z}) / {P}) / ln(2)")
    print(f"   t_cam = -400 * ln({ratio_cam}) / ln(2) = {t_cam:.2f} days.")
    print(f"   Rounded to the nearest day: {a}")

    # --- Part b: Sensor shutdown time ---
    power_sensor_only = x + z
    ratio_sensor = power_sensor_only / P
    # Using the formula: t = -400 * ln(ratio) / ln(2)
    # Since ratio_sensor = 0.05, we calculate -ln(1/0.05) = -ln(20)
    t_sensor = -400 * math.log(ratio_sensor) / math.log(2)
    b = round(t_sensor)
    print(f"\nb. The sensor is turned off when power drops below {power_sensor_only}W.")
    print(f"   The equation is t = -400 * ln(({x} + {z}) / {P}) / ln(2)")
    print(f"   t_sensor = -400 * ln({ratio_sensor}) / ln(2) = {t_sensor:.2f} days.")
    print(f"   Rounded to the nearest day: {b}")

    # --- Part c: Memory usage ---
    # A memory-efficient program would use 2 'frac' variables.
    # size(frac) = size(char) * 3 = 2D * 3 = 6D.
    # Total memory = 2 * 6D = 12D.
    c = 12
    print(f"\nc. The most memory-efficient program would use 2 'frac' variables (6D each).")
    print(f"   Total memory usage = 2 * 6D = {c}D.")

    # --- Part d: Number of exp calls ---
    # The number of calls to exp for ln(v) is ceil(log2(v / 1e-9)).
    # We need to compute ln(2), ln(10), and ln(20).
    calls_ln2 = math.ceil(math.log2(2 / 1e-9))
    calls_ln10 = math.ceil(math.log2(10 / 1e-9))
    calls_ln20 = math.ceil(math.log2(20 / 1e-9))
    d = calls_ln2 + calls_ln10 + calls_ln20
    print(f"\nd. The program calls ln(2), ln(10), and ln(20).")
    print(f"   Number of exp calls for ln(2): {calls_ln2}")
    print(f"   Number of exp calls for ln(10): {calls_ln10}")
    print(f"   Number of exp calls for ln(20): {calls_ln20}")
    print(f"   Total calls to exp = {calls_ln2} + {calls_ln10} + {calls_ln20} = {d}")
    
    # --- Final Answer ---
    print("\nFinal answer in the format a:b:c:d")
    print(f"{a}:{b}:{c}:{d}")


solve_pioneer_problem()
<<<1329:1729:12:100>>>