# Given forward rates
forward_rates = [0.01586388489124628, 0.15527671509596175, 0.10443790375156148, -0.10475743584907427]

# Function to calculate zero rates from forward rates
def calculate_zero_rates(forward_rates):
    zero_rates = [forward_rates[0]]  # The first zero rate is the same as the first forward rate
    for i in range(1, len(forward_rates)):
        # Calculate the zero rate for the i-th year
        product_of_forward_rates = 1.0
        for j in range(i):
            product_of_forward_rates *= (1 + forward_rates[j])
        zero_rate = ((1 + forward_rates[i]) * product_of_forward_rates) ** (1 / (i + 1)) - 1
        zero_rates.append(zero_rate)
    return zero_rates

# Calculate zero rates
zero_rates = calculate_zero_rates(forward_rates)
print(zero_rates)