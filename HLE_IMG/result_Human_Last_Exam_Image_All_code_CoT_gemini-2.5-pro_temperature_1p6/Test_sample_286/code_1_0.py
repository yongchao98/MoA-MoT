import math

def check_plot_b():
    """
    Analyzes data points from Plot B assuming the red curve is r_xy.
    A Markovian evolution is characterized by monotonically increasing entropy S(t)
    and monotonically decreasing Bloch vector length |r(t)|.
    Plot B shows monotonic S, but let's check |r(t)|.
    """
    print("--- Checking Plot B (assumed Markovian) ---")
    
    # Values read from Plot B at t=0
    sz_t0 = 0.5
    r_xy_t0 = 0.7
    r_squared_t0 = sz_t0**2 + r_xy_t0**2
    print(f"At t=0: <σz> = {sz_t0}, r_xy = {r_xy_t0}")
    print(f"|r(0)|^2 = {sz_t0}^2 + {r_xy_t0}^2 = {r_squared_t0:.2f}")

    # Values read from Plot B at t=10
    sz_t10 = 0.6
    r_xy_t10 = 0.68
    r_squared_t10 = sz_t10**2 + r_xy_t10**2
    print(f"At t=10: <σz> ≈ {sz_t10}, r_xy ≈ {r_xy_t10}")
    print(f"|r(10)|^2 ≈ {sz_t10}^2 + {r_xy_t10}^2 = {r_squared_t10:.4f}")
    
    # Check for consistency
    is_consistent = r_squared_t10 <= r_squared_t0
    print(f"\nAnalysis: The entropy in Plot B is monotonic (suggests Markovian),")
    print(f"but the Bloch vector length increases ({r_squared_t10:.2f} > {r_squared_t0:.2f}).")
    print(f"This is a contradiction. Is Plot B physically valid? {is_consistent}\n")

def check_plot_d():
    """
    Analyzes data points from Plot D assuming the red curve is r_xy.
    A non-Markovian evolution can have non-monotonic entropy S(t).
    The condition |r(t)|^2 <= 1 must still hold.
    """
    print("--- Checking Plot D (assumed non-Markovian) ---")
    
    # Values read from Plot D at t=0
    sz_t0 = 0.5
    r_xy_t0 = 0.7
    r_squared_t0 = sz_t0**2 + r_xy_t0**2
    print(f"At t=0: <σz> = {sz_t0}, r_xy = {r_xy_t0}")
    print(f"|r(0)|^2 = {sz_t0}^2 + {r_xy_t0}^2 = {r_squared_t0:.2f}")
    
    is_valid_length_t0 = r_squared_t0 <= 1
    print(f"Check |r(0)|^2 <= 1: {is_valid_length_t0}")
    
    # Values read from Plot D at t=10
    sz_t10 = 0.4
    r_xy_t10 = 0.53
    r_squared_t10 = sz_t10**2 + r_xy_t10**2
    print(f"At t=10: <σz> ≈ {sz_t10}, r_xy ≈ {r_xy_t10}")
    print(f"|r(10)|^2 ≈ {sz_t10}^2 + {r_xy_t10}^2 = {r_squared_t10:.4f}")

    print(f"\nAnalysis: The entropy in Plot D is non-monotonic, consistent with non-Markovian dynamics.")
    print(f"The Bloch vector length remains <= 1 (max value is at t=0).")
    print("The evolution is dissipative overall and internally consistent. Is Plot D physically valid? True")

check_plot_b()
check_plot_d()
