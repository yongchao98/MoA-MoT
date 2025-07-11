import math

# Step 0: Define the initial parameters
engine_rpm = 4500  # Revolutions Per Minute
ecu_mhz = 450       # MegaHertz

# Step 1: Convert engine speed from RPM to RPS (Revolutions Per Second)
seconds_in_minute = 60
engine_rps = engine_rpm / seconds_in_minute

# Step 2: Convert ECU processor speed from MHz to Hz (Cycles Per Second)
mhz_to_hz_multiplier = 1_000_000
ecu_hz = ecu_mhz * mhz_to_hz_multiplier

# Step 3: Calculate the maximum cycles per revolution
# This is the total number of cycles per second divided by the number of revolutions per second.
# (Cycles / Second) / (Revolutions / Second) = Cycles / Revolution
max_cycles_per_revolution = ecu_hz / engine_rps

# --- Output ---
# Display the equation with the specific numbers used in the calculation
print("Equation: (ECU Speed in Hz) / (Engine Speed in RPS)")
print(f"Calculation: ({ecu_mhz} * {mhz_to_hz_multiplier}) / ({engine_rpm} / {seconds_in_minute})")

# Display the final answer
# Using int() to ensure the output is a whole number, as processor cycles are discrete.
print(f"\nThe theoretical maximum number of cycles is: {int(max_cycles_per_revolution)}")