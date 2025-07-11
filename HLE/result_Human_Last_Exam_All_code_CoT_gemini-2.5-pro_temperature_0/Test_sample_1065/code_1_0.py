import math

# Set radius R=1 as the problem asks for raw numbers.
R = 1.0

# Constants
pi = math.pi
sqrt3 = math.sqrt(3)

# The string is composed of two parts:
# 1. An arc on the pumpkin surface from angle 0 to pi/3.
# 2. A vertical segment hanging from the end of the arc.

# Properties of the arc part
# Mass is proportional to length. Let mass density be 1.
m_arc = R * pi / 3
# Center of mass of the arc
# x_cm_arc = (R / (pi/3)) * (1 - cos(pi/3)) = (3*R/pi) * 0.5 = 3*R/(2*pi)
x_cm_arc = 3 * R / (2 * pi)
# y_cm_arc = R + (R / (pi/3)) * sin(pi/3) = R + (3*R/pi) * (sqrt3/2) = R * (1 + 3*sqrt3/(2*pi))
y_cm_arc = R * (1 + 3 * sqrt3 / (2 * pi))

# Properties of the hanging segment part
# Length of the hanging part is R/2
m_seg = R / 2
# The segment hangs from point C(x_c, y_c)
x_c = R * math.sin(pi / 3)
y_c = R + R * math.cos(pi / 3)
# Center of mass of the segment
x_cm_seg = x_c
y_cm_seg = y_c - (m_seg / 2)

# Total mass (proportional to total length)
m_total = m_arc + m_seg

# Calculate the final center of mass coordinates
# X_cm = (m_arc * x_cm_arc + m_seg * x_cm_seg) / m_total
# Y_cm = (m_arc * y_cm_arc + m_seg * y_cm_seg) / m_total

# Horizontal coordinate (X_cm)
x_num = m_arc * x_cm_arc + m_seg * x_cm_seg
x_den = m_total
horizontal_coord = x_num / x_den

# Vertical coordinate (Y_cm)
y_num = m_arc * y_cm_arc + m_seg * y_cm_seg
y_den = m_total
vertical_coord = y_num / y_den

# Print the final coordinates separated by a comma
print(f"{horizontal_coord},{vertical_coord}")
