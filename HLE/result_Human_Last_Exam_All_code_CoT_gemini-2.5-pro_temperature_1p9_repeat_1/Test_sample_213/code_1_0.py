import sys

# This is a helper script to present the answer.
# The shader language is GLSL, not Python.

original_line = 'mat3 Rz=mat3(vec3(cos(a.z),sin(a.z),0),vec3(-sin(a.z),cos(a.z),0),vec3(0,0,1));'

print("The original, unmodified line is:")
# We add the "answer" format tag requested by the user.
print(f"<<<{original_line}>>>")