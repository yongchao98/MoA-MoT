import math

# Given HD-Ratio
hd_ratio = 0.9469858077906315

# Iterate over possible values of TLA and P
for TLA in range(32, 65):  # TLA is typically smaller than P
    for P in range(TLA + 1, 65):  # P should be greater than TLA
        T = math.pow(2, (P - TLA))
        AP = math.pow(T, hd_ratio)
        
        # Check if AP is a reasonable integer
        if AP.is_integer():
            print(f"TLA: {TLA}, P: {P}, AP: {int(AP)}")
            break
    else:
        continue
    break