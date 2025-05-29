import numpy as np

class SalinityPsuCalculatorDouglass2010:
    def __init__(self):
        self.c10 = 42.9

    def from_temperature_c_conductivity_mS_cm(self, temperature_c, conductivity_mS_cm):
        a10 = temperature_c
        b10 = conductivity_mS_cm
        d10 = b10 / self.c10
        g10 = (0.6766097 + 0.0200564 * a10 + 0.0001104259 * a10 ** 2 +
               (-6.9698 * 10 ** -7) * a10 ** 3 + (1.0031 * 10 ** -9) * a10 ** 4)
        e10 = d10 / g10
        f10 = (((a10 - 15) / (1 + 0.0162 * (a10 - 15))) *
               (0.0005 + (-0.0056) * e10 ** 0.5 + (-0.0066) * e10 +
                (-0.0375) * e10 ** 1.5 + (0.0636) * e10 ** 2 +
                (-0.0144) * e10 ** 2.5))
        h10 = (0.008 + (-0.1692 * e10 ** 0.5) + 25.3851 * e10 +
               14.0941 * e10 ** 1.5 + (-7.0261) * e10 ** 2 +
               2.7081 * e10 ** 2.5 + f10)
        return h10

def find_feasible_input(target_salinity):
    calculator = SalinityPsuCalculatorDouglass2010()
    best_diff = float('inf')
    best_temperature = None
    best_conductivity = None

    for temperature in np.arange(0, 35, 0.1):
        for conductivity in np.arange(0, 60, 0.1):
            salinity = calculator.from_temperature_c_conductivity_mS_cm(temperature, conductivity)
            diff = abs(salinity - target_salinity)
            if diff < best_diff:
                best_diff = diff
                best_temperature = temperature
                best_conductivity = conductivity

    return best_temperature, best_conductivity

target_salinity = 21.124100281502077
temperature, conductivity = find_feasible_input(target_salinity)
print({"temperature_c": temperature, "conductivity_mS_cm": conductivity})