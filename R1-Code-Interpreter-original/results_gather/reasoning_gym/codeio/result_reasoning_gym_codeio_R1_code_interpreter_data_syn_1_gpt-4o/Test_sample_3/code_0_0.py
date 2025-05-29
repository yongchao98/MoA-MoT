import math

# Given output
output = [
    {'phase': [-0.7135378710887894, 0.905839606345439], 'magnitude': [6.542064471700567, 11.83242428520853]},
    {'phase': [1.3147148977037437, 1.047913115969739], 'magnitude': [5.269815292499225, 9.907846144313362]},
    {'phase': [2.375467033169534, 2.5726272843237994], 'magnitude': [6.269415378333446, 3.166402434384206]}
]

# Function to calculate real and imaginary parts
def calculate_real_imaginary(phase, magnitude):
    real = magnitude * math.cos(phase)
    imaginary = magnitude * math.sin(phase)
    return real, imaginary

# Calculate real and imaginary parts for each complex image
img_data_list = []
for item in output:
    phase = item['phase']
    magnitude = item['magnitude']
    real_imaginary = []
    for p, m in zip(phase, magnitude):
        real, imaginary = calculate_real_imaginary(p, m)
        real_imaginary.extend([real, imaginary])
    img_data_list.append(real_imaginary)

print(img_data_list)