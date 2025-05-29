import numpy as np
from scipy.optimize import minimize

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

def objective(x):
    ym_per_pix, xm_per_pix = x
    left_curvature, right_curvature = measure_curvature_pixels(ym_per_pix, xm_per_pix)
    target_left = 38.999314499339384
    target_right = 47.443024589235804
    return (left_curvature - target_left)**2 + (right_curvature - target_right)**2

initial_guess = [0.01, 0.01]
result = minimize(objective, initial_guess, bounds=[(0.001, 0.1), (0.001, 0.1)])
ym_per_pix, xm_per_pix = result.x

print({"ym_per_pix": ym_per_pix, "xm_per_pix": xm_per_pix})