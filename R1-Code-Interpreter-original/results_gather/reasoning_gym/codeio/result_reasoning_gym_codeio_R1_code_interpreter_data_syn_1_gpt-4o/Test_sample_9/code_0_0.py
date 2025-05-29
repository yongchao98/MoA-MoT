import math

def find_feasible_inputs(target_hd_ratio):
    for TLA in range(0, 128):
        for P in range(TLA + 1, 129):  # P must be greater than TLA
            T = math.pow(2, (P - TLA))
            AP = math.pow(T, target_hd_ratio)
            if AP.is_integer():
                return TLA, P, int(AP)

target_hd_ratio = 0.9469858077906315
TLA, P, AP = find_feasible_inputs(target_hd_ratio)
print({"TLA": TLA, "P": P, "AP": AP})