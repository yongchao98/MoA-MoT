import numpy as np

def generate_data(ym_per_pix, xm_per_pix):
    np.random.seed(0)
    ploty = np.linspace(0, 719, num=720)
    quadratic_coeff = 3e-4

    leftx = np.array([200 + (y**2)*quadratic_coeff + np.random.randint(-50, high=51) 
                                    for y in ploty])
    rightx = np.array([900 + (y**2)*quadratic_coeff + np.random.randint(-50, high=51) 
                                    for y in ploty])

    leftx = leftx[::-1]
    rightx = rightx[::-1]

    left_fit = np.polyfit(ploty * ym_per_pix, leftx * xm_per_pix, 2)
    right_fit = np.polyfit(ploty * ym_per_pix, rightx * xm_per_pix, 2)
    
    return ploty, left_fit, right_fit

def measure_curvature_pixels(ym_per_pix, xm_per_pix):
    ploty, left_fit, right_fit = generate_data(ym_per_pix, xm_per_pix)

    y_eval = np.max(ploty)
    y_eval *= ym_per_pix

    left_curve_rad = ((1 + (2 * left_fit[0] * y_eval + left_fit[1]) ** 2) ** (3 / 2)) / abs(2 * left_fit[0])
    right_curve_rad = ((1 + (2 * right_fit[0] * y_eval + right_fit[1]) ** 2) ** (3 / 2)) / abs(2 * right_fit[0])

    return left_curve_rad, right_curve_rad

def find_feasible_input():
    target_left_curvature = 38.999314499339384
    target_right_curvature = 47.443024589235804

    # Iterate over a wider range and with smaller increments
    for ym_per_pix in np.arange(0.001, 0.1, 0.001):
        for xm_per_pix in np.arange(0.001, 0.1, 0.001):
            left_curvature, right_curvature = measure_curvature_pixels(ym_per_pix, xm_per_pix)
            if (abs(left_curvature - target_left_curvature) < 0.01 and
                abs(right_curvature - target_right_curvature) < 0.01):
                return {"ym_per_pix": ym_per_pix, "xm_per_pix": xm_per_pix}

    return None

# Execute the function to find the feasible input
feasible_input = find_feasible_input()
print(feasible_input)