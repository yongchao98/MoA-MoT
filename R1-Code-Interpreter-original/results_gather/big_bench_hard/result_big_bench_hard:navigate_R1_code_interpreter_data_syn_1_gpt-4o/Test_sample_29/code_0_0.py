# Calculate net movement
forward_steps = 5
backward_steps = 7 + 10
net_movement = forward_steps - backward_steps

# Check if we return to the starting point
return_to_start = (net_movement == 0)

print("Yes" if return_to_start else "No")