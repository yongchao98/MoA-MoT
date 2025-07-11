# Based on the principles of TIG (GTAW) welding for Inconel 718 repair:
# 1. Voltage Requirement: The specified arc gap of 6 mm is large for TIG welding and
#    requires a correspondingly high voltage to sustain a stable arc.
#    A voltage of approximately 15 V is necessary for this arc length.
# 2. Current Requirement: The task is a material build-up on a turbine blade tip.
#    This requires sufficient energy to melt both the base metal (Inconel 718)
#    and the filler material to ensure proper fusion. A current of 25 A is too low for this
#    application and would not create a stable weld pool. A current of 100 A is a
#    standard and appropriate level of energy for this type of repair.

# Therefore, the most suitable parameters are selected.
selected_current = 100 # Amps
selected_voltage = 15  # Volts
travel_speed = 0.5     # mm/s
arc_gap = 6            # mm

print("Welding Procedure Specification for Root Pass:")
print("---------------------------------------------")
print(f"Material: Inconel 718")
print(f"Process: Manual TIG (GTAW) Repair")
print(f"Arc Gap: {arc_gap} mm")
print(f"Travel Speed: {travel_speed} mm/s")
print("\nRecommended Parameters for a stable repair:")
print(f"Current: {selected_current} A")
print(f"Voltage: {selected_voltage} V")
