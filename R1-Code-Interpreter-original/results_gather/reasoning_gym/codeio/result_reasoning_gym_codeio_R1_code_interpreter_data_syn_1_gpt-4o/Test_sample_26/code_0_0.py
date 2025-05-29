# Constants
D_LEN = 7  # Diatonic length
C_LEN = 12  # Chromatic length

# Input values
x = [-50, -48, 89]
y = [64, -49]

# Compute the tonal sum
d_sum = x[0] + y[0]
c_sum = x[1] + y[1]

# Normalize the diatonic and chromatic values
d_val = d_sum % D_LEN
d_oct = d_sum // D_LEN
c_val = c_sum % C_LEN

# Calculate the resulting octave
o_val = x[2] + d_oct

# Resulting tonal value
result = [d_val, c_val, o_val]

# Output the result
print(result)