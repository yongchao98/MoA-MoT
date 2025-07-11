import sys

# Define the values from the circuit diagram
V = 41.0  # Voltage of the source in Volts
R11 = 11.0 # Resistance in Ohms
R29 = 29.0 # Resistance in Ohms

# --- Circuit Analysis ---
# The wire carrying current I connects the 41V source directly to the node
# where the 11 Ohm and 29 Ohm resistors meet.
# This means the 11 Ohm and 29 Ohm resistors are in parallel, and the
# voltage across them is 41V.
# The current I is the total current flowing into this parallel combination.
# By Kirchhoff's Current Law, I is the sum of the currents through R11 and R29.
# The other resistors (76, 8, 14 Ohms) do not affect the value of I in this configuration.

# Calculate the current through each relevant resistor
I_11 = V / R11
I_29 = V / R29

# Calculate the total current I
I_total = I_11 + I_29

# --- Output the Calculation ---
print("To find the value of the current I, we can use Ohm's law.")
print("The current I is the sum of the currents through the 11 Ohm and 29 Ohm resistors, which are in parallel across the 41V source.")
print("\nThe equation for the total current I is:")
print("I = (Voltage / Resistance_1) + (Voltage / Resistance_2)")
print("\nSubstituting the values from the circuit diagram:")
# Using sys.stdout.write to avoid the default space of print
sys.stdout.write("I = " + str(V) + "V / " + str(R11) + " Ohms + " + str(V) + "V / " + str(R29) + " Ohms\n")
sys.stdout.write("I = " + str(I_11) + " A + " + str(I_29) + " A\n")
sys.stdout.write("I = " + str(I_total) + " A\n")