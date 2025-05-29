solution = [5, 15, 20, 15, 10, 10, 36, 50, 45, 96, 3, 24, 16]
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M']

# Verify all conditions
F, G = solution[5], solution[6]
J, M = solution[9], solution[12]
H, L = solution[7], solution[11]
B = solution[1]
A, C = solution[0], solution[2]
I = solution[8]
D, E = solution[3], solution[4]

print(f"Verification:")
print(f"1. F + G = {F + G} (should be 46)")
print(f"2. J - M = {J - M} (should be 80)")
print(f"3. J - F = {J - F} (should be 86)")
print(f"4. H - L = {H - L} (should be 26)")
print(f"5. B + M = {B + M} (should be 31)")
print(f"6. B = {B}, 1.5F = {1.5*F}")
print(f"7. C = {C}, 4.0A = {4.0*A}")
print(f"8. L = {L}, 1.5M = {1.5*M}")
print(f"9. G = {G}, 3.6F = {3.6*F}")
print(f"10. D = {D}, 1.5E = {1.5*E}")
print(f"11. L - I = {L - I} (should be -21)")