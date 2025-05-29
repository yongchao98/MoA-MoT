import math
import random

class BoxMueller:
    call_count = 0
    second_random_variable = 0.0
    
    @classmethod
    def randn(cls):
        cls.call_count += 1
        if cls.call_count % 2 == 0:
            return cls.second_random_variable
        
        u1 = 1.0 - random.random()
        u2 = 1.0 - random.random()
        c = math.sqrt(-2.0 * math.log(u1))
        cls.second_random_variable = c * math.sin(2 * math.pi * u2)
        return c * math.cos(2 * math.pi * u2)

def main_solution(seed, freq_1, freq_2):
    random.seed(seed)
    
    white_noise = [BoxMueller.randn() for _ in range(16)]
    
    inner_product_1 = sum(math.cos(freq_1 * 2.0 * math.pi * (i + 0.5) / 16.0) * white_noise[i] for i in range(16))
    inner_product_2 = sum(math.cos(freq_2 * 2.0 * math.pi * (i + 0.5) / 16.0) * white_noise[i] for i in range(16))
    
    return {
        "inner_product_1": inner_product_1,
        "inner_product_2": inner_product_2
    }

# Given input
input_data = {'seed': 914215, 'freq_1': 3, 'freq_2': 6}

# Calculate the output
output = main_solution(input_data['seed'], input_data['freq_1'], input_data['freq_2'])

# Print the output
print(output)