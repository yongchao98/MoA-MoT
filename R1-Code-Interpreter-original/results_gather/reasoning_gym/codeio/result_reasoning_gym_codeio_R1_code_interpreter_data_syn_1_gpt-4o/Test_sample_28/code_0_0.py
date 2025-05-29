import numpy as np

def function_one(x, sigma):
    return (x**2) * np.exp((-1/2)*((x/sigma)**2))

def function_one_deriv(x, sigma):
    return (((x**3)*np.exp(-(x**2)/(2 * (sigma**2))))/(2 * sigma))+(2*x*np.exp(-(x**2)/(2 * (sigma**2))))

def function_two(x, sigma):
    return np.exp((-1/2)*((x/sigma)**2))

def function_two_deriv(x, sigma):
    return (np.exp(-(x**2)/(2 * (sigma**2))))/(2 * sigma)

def integration_F1(xi, xf, n, sigma):
    dx = (xf - xi)/n
    total1 = 0
    for i in range(n):
        x = xi + i * dx
        f = function_one(x, sigma)
        area = f * dx
        total1 += area
    return total1

def integration_F2(xi, xf, n, sigma):
    dx = (xf - xi)/n
    total2 = 0
    for i in range(n):
        x = xi + i * dx
        f = function_two(x, sigma)
        area = f * dx
        total2 += area
    return total2

def uncertainty_F1(xi, xf, n, sigma):
    x = np.linspace(xi, xf, n)
    dxdt_1 = function_one_deriv(x, sigma)
    M1 = dxdt_1.max()
    return (1/2) * M1 * (((xf - xi)**2)/n)

def uncertainty_F2(xi, xf, n, sigma):
    x = np.linspace(xi, xf, n)
    dxdt_2 = function_two_deriv(x, sigma)
    M2 = dxdt_2.max()
    return (1/2) * M2 * (((xf - xi)**2)/n)

def inertia(M, xi, xf, n, sigma):
    return M * (integration_F1(xi, xf, n, sigma)/integration_F2(xi, xf, n, sigma))

def uncert_inertia(M, xi, xf, n, sigma):
    return M * ((uncertainty_F1(xi, xf, n, sigma)/integration_F1(xi, xf, n, sigma))+(uncertainty_F2(xi, xf, n, sigma)/integration_F2(xi, xf, n, sigma)))* inertia(M, xi, xf, n, sigma)

def main_solution(M, L, sigma, n):
    xi = -L/2
    xf = L/2
    inertia_value = inertia(M, xi, xf, n, sigma)
    uncertainty_value = uncert_inertia(M, xi, xf, n, sigma)
    return {"inertia": inertia_value, "uncertainty": uncertainty_value}

# Given input
input_data = {'M': 2.9058664904667655, 'L': 1.193615739757427, 'sigma': 0.5363370748982019, 'n': 595}
result = main_solution(**input_data)
print(result)